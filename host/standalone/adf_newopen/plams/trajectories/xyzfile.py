#!/usr/bin/env python

import numpy
from ..mol.molecule import Molecule
from ..tools.geometry import cell_shape
from .trajectoryfile import TrajectoryFile

__all__ = ['XYZTrajectoryFile','create_xyz_string']

class XYZTrajectoryFile (TrajectoryFile) :
        """
        Class representing an XYZ file containing a molecular trajectory

        An instance of this class has the following attributes:

        *   ``file_object`` -- A Python :py:class:`file` object, referring to the actual XYZ file
        *   ``position``    -- The frame to which the cursor is currently pointing in the XYZ file
        *   ``mode``        -- Designates whether the file is in read or write mode ('r' or 'w')
        *   ``ntap``        -- The number of atoms in the molecular system (needs to be constant throughout)
        *   ``elements``    -- The elements of the atoms in the system (needs to be constant throughout)

        An |XYZTrajectoryFile| object behaves very similar to a regular file object.
        It has read and write methods (:meth:`read_next` and :meth:`write_next`) 
        that read and write from/to the position of the cursor in the ``file_object`` attribute. 
        If the file is in read mode, an additional method :meth:`read_frame` can be used that moves
        the cursor to any frame in the file and reads from there.
        The amount of information stored in memory is kept to a minimum, as only information from the current frame
        is ever stored.

        Reading and writing to and from the files can be done as follows::

            >>> from scm.plams import XYZTrajectoryFile

            >>> xyz = XYZTrajectoryFile('old.xyz')
            >>> mol = xyz.get_plamsmol()

            >>> xyzout = XYZTrajectoryFile('new.xyz',mode='w')

            >>> for i in range(xyz.get_length()) :
            >>>     crd,cell = xyz.read_frame(i,molecule=mol)
            >>>     xyzout.write_next(molecule=mol)

        The above script reads information from the XYZ file ``old.xyz`` into the |Molecule| object ``mol``
        in a step-by-step manner.
        The |Molecule| object is then passed to the :meth:`write_next` method of the new |XYZTrajectoryFile|
        object corresponding to the new xyz file ``new.xyz``.

        The exact same result can also be achieved by iterating over the instance as a callable

            >>> xyz = XYZTrajectoryFile('old.xyz')
            >>> mol = xyz.get_plamsmol()

            >>> xyzout = XYZTrajectoryFile('new.xyz',mode='w')

            >>> for crd,cell in xyz(mol) :
            >>>     xyzout.write_next(molecule=mol)

        This procedure requires all coordinate information to be passed to and from the |Molecule| object
        for each frame, which can be time-consuming.
        It is therefore also possible to bypass the |Molecule| object when reading through the frames::

            >>> xyz = XYZTrajectoryFile('old.xyz')

            >>> xyzout = XYZTrajectoryFile('new.xyz',mode='w')
            >>> xyzout.set_elements(xyz.get_elements())

            >>> for crd,cell in xyz :
            >>>     xyzout.write_next(coords=crd)
            >>> xyzout.close()

        By default the write mode will create a minimal version of the XYZ file, containing only elements
        and coordinates. 
        Additional information can be written to the file by supplying additional arguments
        to the :meth:`write_next` method. 
        The additional keywords `step` and `energy` trigger the writing of a remark containing
        the molecule name, the step number, the energy, and the lattice vectors.

            >>> mol = Molecule('singleframe.xyz')

            >>> xyzout = XYZTrajectoryFile('new.xyz',mode='w')
            >>> xyzout.set_name('MyMol')

            >>> xyzout.write_next(molecule=mol, step=0, energy=5.)
        """
        def __init__ (self, filename, mode='r', fileobject=None, ntap=None) :
                """
                Initiates an XYZTrajectoryFile object

                * ``filename``   -- The path to the XYZ file
                * ``mode``       -- The mode in which to open the XYZ file ('r' or 'w')
                * ``fileobject`` -- Optionally, a file object can be passed instead (filename needs to be set to None)
                * ``ntap``       -- If the file is in write mode, the number of atoms needs to be passed here
                """
                TrajectoryFile.__init__(self, filename, mode, fileobject, ntap)

                # XYZ specific attributes
                self.name = 'PlamsMol'

                # Required setup before frames can be read/written
                if self.mode == 'r' :
                        self._read_header()

        def set_name (self, name) :
                """
                Sets the name of the system, in case an extensive write is requested

                *   ``name`` -- A string containing the name of the molecule
                """
                self.name = name

        def _read_header (self) :
                """
                Set up info required for reading frames
                """
                line = self.file_object.readline()
                self.ntap = int(line.split()[0])
                if self.coords.shape == (0,3) :
                        self.coords = numpy.zeros((self.ntap,3)) 

                self.file_object.readline()

                elements = []
                for i in range(self.ntap) :
                        line = self.file_object.readline()
                        elements.append(line.split()[0])
                self.elements = elements

                self.file_object.seek(0)

        def read_next (self, molecule=None, read=True) :
                """
                Reads coordinates from the current position of the cursor and returns it

                * ``molecule`` -- |Molecule| object in which the new coordinates need to be stored
                * ``read``     -- If set to False the cursor will move to the next frame without reading
                """
                if not read and not self.firsttime :
                        return self._move_cursor_without_reading()

                cell = None
                # Read the coordinates
                crd = self._read_coordinates(molecule)
                if crd is None :
                        return None, None       # End of file is reached

                if self.firsttime :
                        self.firsttime = False

                self.position += 1
                
                return self.coords, cell

        def _read_coordinates (self, molecule) :
                """
                Read the coordinates from file, and place them in the molecule
                """
                cell = None
                for i in range(2) :
                        line = self.file_object.readline()
                        if len(line) == 0 :
                                return None           # End of file is reached
                for i in range(self.ntap) :
                        line = self.file_object.readline()
                        self.coords[i,:] = [float(w) for w in line.split()[1:4]]

                if isinstance(molecule,Molecule) :
                        self._set_plamsmol(self.coords, cell, molecule)

                return self.coords

        def _is_endoffile (self) :
                """
                If the end of file is reached, return coords and cell as None
                """
                end = False
                for i in range(self.ntap+2) :
                        line = self.file_object.readline()
                        if len(line) == 0 :
                                end = True
                                break
                return end

        def write_next (self,coords=None,molecule=None,cell=[0.,0.,0.],energy=None,step=None,conect=None) :
                """
                Write frame to next position in trajectory file

                * ``coords``   -- A list or numpy array of (``ntap``,3) containing the system coordinates
                * ``molecule`` -- A molecule object to read the molecular data from
                * ``cell``     -- A set of lattice vectors or cell diameters
                * ``energy``   -- An energy value to be written to the remark line
                * ``step``     -- A frame number for the current frame
                * ``conect``   -- A dictionary containing connectivity info (not used)

                .. note::

                        Either ``coords`` or ``molecule`` are mandatory arguments
                """
                if isinstance(molecule,Molecule) :
                        coords, cell, elements = self._read_plamsmol(molecule)[:3]
                        if self.position == 0 :
                                self.elements = elements
                cell = self._convert_cell(cell)

                self._write_moldata(coords, cell, energy, step)

                self.position += 1

        def _write_moldata (self, coords, cell, energy, step) :
                """
                Write all molecular info to file
                """
                if step is not None or energy is not None :
                        if step is None :
                                step = self.position
                        if energy is None :
                                energy = 0.
                        box = None
                        if cell is not None :
                                #box = PDBMolecule().box_from_vectors(cell)
                                box = cell_shape(cell)
                        block = create_xyz_string(self.elements,coords,energy,box,step,self.name)
                else :
                        block = create_xyz_string(self.elements, coords)
                self.file_object.write(block)

        def _rewind_to_first_frame(self) :
                """
                Rewind the file to the first frame
                """
                self.file_object.seek(0)
                self.firsttime = True
                self.position = 0

        def _rewind_n_frames(self,nframes) :
                """
                Rewind the file by nframes frames
                """
                new_frame = self.position - nframes
                self._rewind_to_first_frame()
                for i in range(new_frame) :
                        self.read_next(read=False)


def create_xyz_string (elements, coords, energy=None, box=None, step=None, name='PlamsMol') :
        """
        Write an XYZ file based on the elements and the coordinates of the atoms
        """
        block = '%8i\n'%(len(elements))
        if step is not None :
                if energy is None : energy = 0.
                comment = '%-40s%6i %16.6f'%(name, step, energy)
                if box is not None :
                        for value in box :
                                comment += '%7.2f'%(value)
                block += comment
        block += '\n'
        for el,crd in zip(elements,coords) :
                block += '%8s '%(el)
                for x in crd :
                        block += '%20.10f '%(x)
                block += '\n'
        return block
