from .scmjob import SCMJob, SCMResults
from ...core.errors import ResultsError
from ...core.functions import log
from ...tools.units import Units

__all__ = ['DFTBJob', 'DFTBResults']



class DFTBResults(SCMResults):
    _kfext = '.rkf'
    _rename_map = {'dftb.rkf':'$JN'+_kfext}


    def get_main_molecule(self):
        """get_main_molecule()
        Return a |Molecule| instance based on the ``Molecule`` section in the main KF file (``$JN.rkf``).

        For runs with multiple geometries (geometry optimization, transition state search, molecular dynamics) this is the **final** geometry. In such a case, to access the initial (or any intermediate) coordinates please extract them from section ``History``, variable ``Coords(i)`` or use :meth:`get_input_molecule`. Mind the fact that all coordinates written by DFTB to ``.rkf`` file are in bohr::

            mol = results.get_molecule(section='History', variable='Coords(1)', unit='bohr')
        """
        ret = self.get_molecule(section='Molecule', variable='Coords', unit='bohr')
        ret.properties.charge = self.readkf('Molecule', 'Charge')
        if ('Molecule', 'LatticeVectors') in self._kf:
            lattice = self.readkf('Molecule', 'LatticeVectors')
            lattice = Units.convert(lattice, 'bohr', 'angstrom')
            ret.lattice = [tuple(lattice[i:i+3]) for i in range(0,len(lattice),3)]
        return ret


    def get_input_molecule(self):
        """get_input_molecule()
        Return a |Molecule| instance with initial coordinates.

        All data used by this method is taken from ``$JN.rkf`` file. The ``molecule`` attribute of the corresponding job is ignored.
        """
        if ('History', 'nEntries') in self._kf:
            return self.get_molecule(section='History', variable='Coords(1)', unit='bohr')
        return self.get_main_molecule()


    def get_history_molecules(self):
        """get_history_molecules()
        Return a list of |Molecule| instances with all Structures from the ``History`` section of the ``.rkf`` file.
        If the structures have lattice information, it is written into the corresponding |Molecule| instance.

        All data used by this method is taken from ``$JN.rkf`` file. If no ``History`` section is available an empty list is returned and a level 5 log entry appears.
        """
        if not ('History', 'nEntries') in self._kf:
            log('No History section found in rkf during get_history_molecules, returning an empty list.',level=5)
            return []

        nEntries = self.readkf('History', 'nEntries')
        #list of molecules to return
        history = []
        #iterate over entries
        for i in range(1,nEntries+1):
            mol = self.get_molecule('History','Coords('+str(i)+')')
            #read lattice if there is any
            if ('History','LatticeVectors('+str(i)+')') in self._kf:
                #LatticeVectors has to be length dimensions*3, so we can do this
                lattice = Units.convert(self.readkf('History','LatticeVectors('+str(i)+')'), 'bohr', 'angstrom')
                mol.lattice = [tuple(lattice[j:j+3]) for j in range(0,len(lattice),3)]
            history.append(mol)
        return history


    def get_energy(self, unit='au'):
        """get_energy(unit='au')
        Return DFTB final energy, expressed in *unit*.
        """
        prop = self.get_properties()
        if 'DFTB Final Energy' in prop:
            return Units.convert(prop['DFTB Final Energy'], 'au', unit)
        raise ResultsError("'DFTB Final Energy' not present in 'Properties' section of {}".format(self._kfpath()))


    def get_dipole_vector(self, unit='au'):
        """get_dipole_vector(unit='au')
        Return the dipole vector, expressed in *unit*.
        """
        prop = self.get_properties()
        if 'Electric Dipole Moment' in prop:
            return Units.convert(prop['Electric Dipole Moment'], 'au', unit)
        raise ResultsError("'Electric Dipole Moment' not present in 'Properties' section of {}".format(self._kfpath()))


    def get_gradients(self):
        """get_gradients()
        Return the list of atomic gradients, expressed in atomic units. Returned value is a list of vectors (lists of lenght 3).
        """
        prop = self.get_properties()
        if 'Generic Gradient' in prop:
            grad = prop['Generic Gradient']
            return [grad[i:i+3] for i in range(0,len(grad),3)]
        raise ResultsError("'Generic Gradient' not present in 'Properties' section of {}".format(self._kfpath()))


    def _int2inp(self):
        """_int2inp()
        In DFTB the internal order is always the same as the input order. Return an identity permutation of length equal to the number of atoms.
        """
        return list(range(1, 1+self.readkf('Molecule', 'nAtoms')))


    def _atomic_numbers_input_order(self):
        """_atomic_numbers_input_order()
        Return a list of atomic numbers, in the input order.
        """
        return self.readkf('Molecule', 'AtomicNumbers')



class DFTBJob(SCMJob):
    _result_type = DFTBResults
    _command = 'dftb'
    _top = ['units', 'task']
    _subblock_end = 'end'


    def _serialize_mol(self):
        s = self.settings.input

        if len(self.molecule.lattice) in [1,2] and self.molecule.align_lattice():
            log("The lattice supplied for job {} did not follow the convention required by DFTB. I rotated the whole system for you. You're welcome".format(self._full_name()), 3)

        for i,atom in enumerate(self.molecule):
            s.system.atoms['_'+str(i+1)] = atom.str(symbol=self._atom_symbol(atom), space=18, decimal=10)

        if self.molecule.lattice:
            for i,vec in enumerate(self.molecule.lattice):
                s.system.lattice['_'+str(i+1)] = '{:16.10f} {:16.10f} {:16.10f}'.format(*vec)


    def _remove_mol(self):
        s = self.settings.input
        if 'system' in s:
            if 'atoms' in s[system]:
                del s.system.atoms
            if 'lattice' in s.system:
                del s.system.lattice
