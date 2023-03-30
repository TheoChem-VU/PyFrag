#!/usr/bin/env python

import os
import numpy
from..mol.molecule import Bond
from ..core.errors import PlamsError

__all__ = ['TrajectoryFile']

class TrajectoryFile (object) :
        """
        Abstract class that represents a generic trajectory file
        """
        def __init__ (self, filename, mode='r', fileobject=None, ntap=None) :
                """
                Would create a generic trajectory  file object

                * ``filename``   -- The path to the DCD file
                * ``mode``       -- The mode in which to open the DCD file ('rb' or 'wb')
                * ``fileobject`` -- Optionally, a file object can be passed instead (filename needs to be set to None)
                * ``ntap``       -- If the file is in write mode, the number of atoms needs to be passed here
                """
                self.position = 0
                if filename is not None :
                        fileobject = open(filename,mode)
                self.file_object = fileobject
                if self.file_object is not None :
                        self.mode = self.file_object.mode

                self.ntap = 0
                if ntap is not None :
                        self.ntap = ntap
                self.firsttime = True
                self.coords = numpy.zeros((self.ntap,3))                # Only for reading purposes, to avoid creating the array each time

                # PLAMS molecule related settings
                self.elements = ['H']*self.ntap
                self.current_molecule = None
                self.store_molecule = True # Even if True, the molecule attribute is only stored during iteration

        def __iter__ (self) :
                """
                Magic method that makes this an iterator

                Note: Iterates starting from current cursor position
                """
                if self.mode[0] == 'r' :
                        if self.current_molecule is None and self.store_molecule :
                                self.current_molecule = self.get_plamsmol()
                        for i in range(len(self)-self.position) :
                                crd,cell = self.read_next(molecule=self.current_molecule)
                                yield crd,cell
                        self.current_molecule = None
                else :
                        raise PlamsError('Cannot iterate over writable trajectoryfile')

        @property
        def molecule (self) :
                """
                Returns the current molecule, which exists only when the object is used as an iterator
                """
                if not self.store_molecule :
                        raise PlamsError ('Try setting the store_molecule attribute to True')
                elif self.current_molecule is None :
                        raise PlamsError ('molecule attribute only exists during iteration')
                return self.current_molecule

        def __call__ (self, molecule=None, read=True) :
                """
                Magic method that makes an instance of this class into a callable
                """
                if self.mode[0] == 'r' :
                        for i in range(len(self)-self.position) :
                                crd,cell = self.read_next(molecule, read)
                                yield crd,cell
                else :
                        raise PlamsError('Cannot iterate over writable trajectoryfile')

        def __len__ (self) :
                """
                Magic method that allows checking the length of an iterator
                """
                return self.get_length()

        def __enter__ (self) :
                """
                Magic method that allows use with *with* statement
                """
                return self

        def __exit__ (self, exc_type, exc_val, exc_tb) :
                """
                Magic method that allows use with *with* statement
                """
                self.close()

        def __del__ (self) :
                """
                If self.file_object is a regular fileobject, close it
                """
                try :
                        self.file.close()
                except AttributeError :
                        pass

        def get_elements (self) :
                """
                Get the elements attribute
                """
                return self.elements

        def set_elements (self, elements) :
                """
                Sets the elements attribute (needed in write mode).

                *   ``elements`` -- A list containing the element symbol of each atom
                """
                self.elements = elements

        def _read_header (self) :
                """
                Set up info required for reading frames

                Sets self.ntap and self.coords to proper value/size
                """
                pass

        def get_plamsmol (self) :
                """
                Extracts a PLAMS molecule object from the XYZ trajectory file
                """
                from scm.plams import Molecule
                oldposition = self.position
                self.rewind()
                coords, cell = self.read_next()
                plamsmol = Molecule.from_elements(self.elements)
                plamsmol.from_array(coords)

                # Return to original position
                self.rewind()
                for i in range(oldposition) :
                        self.read_next(read=False)

                return plamsmol

        def close (self) :
                """
                Cleanly close and garbage collect the file
                """
                del(self)

        def read_next (self, molecule=None, read=True) :
                """
                Reads the relevant info from current frame
                """
                if not read and not self.firsttime :
                        return self._move_cursor_without_reading()

                # Read the coordinates and cell info from self.file_object
                self.coords = None
                self.cell = None

                # Finalize
                if self.firsttime :
                        self.firsttime = False

                # Place the values into the provided molecule object
                if isinstance(molecule,Molecule) :
                        self._set_plamsmol(self.coords,cell,molecule)

                self.position += 1
                return self.coords, cell

        def _move_cursor_without_reading (self) :
                """
                Move the cursor forward and return coords as empty array
                """
                cell = numpy.zeros((3,3))
                if self._is_endoffile() :
                        return None, None
                self.position += 1
                coords = numpy.array([])
                return coords, cell

        def _is_endoffile (self) :
                """
                Reads and checks If the end of file is reached.
                """
                pass

        def read_frame (self, frame, molecule=None) :
                """
                Reads the relevant info from frame ``frame`` and returns it, or stores it in ``molecule``

                * ``frame``    -- The frame number to be read from the file
                * ``molecule`` -- |Molecule| object in which the new coordinates need to be stored
                """
                if frame < self.position :
                        nframes = abs(frame - self.position)
                        self.rewind(nframes)
                steps = frame - self.position
                
                for i in range(steps) :
                        crd,cell = self.read_next(read=False)
                        if crd is None :
                                break
                
                crd,cell = self.read_next(molecule)
                if crd is None :
                        print('Not enough frames!')
                
                return crd,cell

        def _set_plamsmol (self, coords, cell, plamsmol, bonds=None) :
                """
                If molecule objects have been passed as arguments, update their coordinates and lattice
                """
                plamsmol.from_array(coords)
                if cell is not None :
                        if cell[0,0] > 0. :
                                plamsmol.lattice = cell.tolist()
                if bonds is not None :
                        plamsmol.delete_all_bonds()
                        for bond in bonds :
                                b = Bond(plamsmol[bond[0]],plamsmol[bond[1]])
                                plamsmol.add_bond(b)

        def write_next (self, coords=None, molecule=None, cell=[0.,0.,0.], energy=0., step=None, conect=None) :
                """
                Write the coordinates to the next available spot in the file
                """
                if isinstance(molecule,Molecule) :
                        coords, cell = self._read_plamsmol(molecule)[:2]
                cell = self._convert_cell(cell)
                # Write to self.file_object
                self.position += 1

        def _convert_cell (self, cell) :
                """
                Check format of cell, and convert to vectors if needed
                """
                if cell is None :
                        return cell
                if len(cell)==3 and isinstance(cell[0],float) or isinstance(cell[0],numpy.float64) :
                        cell = numpy.diag(cell)
                else :
                        cell = numpy.array(cell)
                if cell[0,0] == 0. :
                        cell = None
                return cell

        def _read_plamsmol (self, plamsmol) :
                """
                Read the coordinates and cell vectors from the molecule objects, if provided
                """
                coords = plamsmol.as_array()
                cell = [0.,0.,0.]
                if len(plamsmol.lattice) > 0 : cell = plamsmol.lattice
                elements = [at.symbol for at in plamsmol.atoms]
                # Get the connection table
                plamsmol.set_atoms_id()
                conect = {}
                for iat,atom in enumerate(plamsmol.atoms) :
                        neighbors = plamsmol.neighbors(atom)
                        if len(neighbors) > 0 :
                                conect[iat+1] = [neighbor.id for neighbor in neighbors]
                # Get the properties
                props = [at.properties.suffix if 'suffix' in at.properties else '' for at in plamsmol.atoms]
                if sum([len(s) for s in props]) == 0 : props = None
                return coords, cell, elements, conect, props

        def rewind (self,nframes=None) :
                """ 
                Rewind the file either by ``nframes`` or to the first frame

                *   ``nframes`` -- The number of frames to rewind
                """
                if nframes is None or nframes == self.position :
                        # Go back to the beginning
                        self._rewind_to_first_frame()
                elif nframes > self.position :
                        raise PlamsError('Trying to rewind too much!')
                else :  
                        # Go back nframes geometries
                        self._rewind_n_frames(nframes)

        def _rewind_to_first_frame(self) :
                """
                Rewind the file to the first frame
                """
                self.file_object.seek(0)
                self.firsttime = True
                self.position = 0
                self._read_header()
                
        def _rewind_n_frames(self,nframes) :
                """
                Rewind the file by nframes frames
                """
                pass

        def get_length (self) :
                """
                Get the number of frames in the file
                """
                oldposition = self.position
                while 1 :
                        crd,cell = self.read_next(read=False)
                        if crd is None :
                                break
                #size = self.position-1
                size = self.position

                # Go back to the original position
                self.rewind()
                for i in range(oldposition) :
                        crd,cell = self.read_next(read=False)

                return size

        def read_last_frame (self, molecule=None) :
                """
                Reads the last frame from the file
                """
                step = 0
                while 1 :
                        crd,cell = self.read_next(read=False)
                        if crd is None :
                                break
                        step += 1

                # If the current position was already the end of the file, this fails, so a failsafe is built in:                
                if step == 0 :
                        step = self.position

                self.rewind(1)
                crd,cell = self.read_next(molecule=molecule)
                return crd,cell
