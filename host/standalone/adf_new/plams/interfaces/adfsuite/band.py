from .scmjob import SCMJob, SCMResults
from ...core.errors import ResultsError
from ...core.functions import log
from ...tools.units import Units

__all__ = ['BANDJob', 'BANDResults']


class BANDResults(SCMResults):
    _kfext = '.runkf'
    _rename_map = {'RUNKF':'$JN'+_kfext}


    def get_properties(self):
        """get_properties()
        Return a dictionary with all the entries from ``Properties`` section in the main KF file (``$JN.rkf``).
        """
        n = self.readkf('Properties', 'nEntries')
        ret = {}
        for i in range(1, n+1):
            tp = self.readkf('Properties', 'Type({})'.format(i)).strip()
            stp = self.readkf('Properties', 'Subtype({})'.format(i)).strip()
            val = self.readkf('Properties', 'Value({})'.format(i))
            if stp.endswith(' (a.u.)'):
                stp = stp[:-7]
            key = '{} {}'.format(stp, tp.lower())
            ret[key] = val
        return ret


    def get_main_molecule(self):
        """get_main_molecule()
        Return a |Molecule| instance based on the ``Molecule`` section in the main KF file (``$JN.runkf``).

        For runs with multiple geometries (geometry optimization, transition state search, molecular dynamics) this is the **final** geometry. In such a case, to access the initial (or any intermediate) coordinates please use :meth:`get_input_molecule` or extract coordinates from section ``History``, variables ``xyz0``, ``xyz1`` and so on. Mind the fact that all coordinates written by BAND to ``History`` section are in bohr and internal atom order::

            mol = results.get_molecule(section='History', variable='xyz0', unit='bohr', internal=True)
        """
        ret = self.get_molecule(section='Molecule', variable='Coords', unit='bohr', internal=True)
        ret.properties.charge = self.readkf('Molecule', 'Charge')
        if ('Molecule', 'LatticeVectors') in self._kf:
            lattice = self.readkf('Molecule', 'LatticeVectors')
            lattice = Units.convert(lattice, 'bohr', 'angstrom')
            ret.lattice = [tuple(lattice[i:i+3]) for i in range(0,len(lattice),3)]
        return ret


    def get_input_molecule(self):
        """get_input_molecule()
        Return a |Molecule| instance with initial coordinates.

        All data used by this method is taken from ``$JN.runkf`` file. The ``molecule`` attribute of the corresponding job is ignored.
        """
        if ('History', 'nr of geometries') in self._kf:
            ret = self.get_molecule(section='History', variable='xyz0', unit='bohr', internal=False)
            lattice = self.readkf('History', 'lattice0')
            lattice = Units.convert(lattice, 'bohr', 'angstrom')
            lattice = [lattice[i:i+3] for i in range(0,len(lattice),3)]
            while len(lattice) > 0 and not any(lattice[-1]):
                lattice.pop()
            ret.lattice = [tuple(i) for i in lattice]
            return ret
        return self.get_main_molecule()


    def get_energy(self, unit='au'):
        """get_energy(unit='au')
        Return final bond energy, expressed in *unit*.
        """
        return self._get_single_value('Bond energies', 'final bond energy', unit)


    def get_fermi_energy(self, unit='au'):
        """get_fermi_energy(unit='au')
        Return Fermi energy, expressed in *unit*.
        """
        return self._get_single_value('BandStructure', 'FermiEnergy', unit)


    def get_band_gap(self, unit='au'):
        """get_band_gap(unit='au')
        Return the value of the band gap, expressed in *unit*.
        """
        return self._get_single_value('BandStructure', 'BandGap', unit)


    def get_dipole_vector(self, unit='au'):
        """get_dipole_vector(unit='au')
        Return the dipole vector, expressed in *unit*.
        """
        prop = self.get_properties()
        if 'Dipole moment vector' in prop:
            return Units.convert(prop['Dipole moment vector'], 'au', unit)
        raise ResultsError("'Dipole moment vector' not present in 'Properties' section of {}".format(self._kfpath()))


    def _int2inp(self):
        """_int2inp()
        Get mapping from the internal atom order to the input atom order.
        """
        return self.readkf('geometry', 'Atom map new order')


    def _atomic_numbers_input_order(self):
        """_atomic_numbers_input_order()
        Return a list of atomic numbers, in the input order.
        """
        mapping = self._int2inp()
        atnums = self.readkf('Molecule', 'AtomicNumbers')
        return [atnums[mapping[i]-1] for i in range(len(atnums))]



class BANDJob(SCMJob):
    _result_type = BANDResults
    _command = 'band'


    def _serialize_mol(self):
        s = self.settings.input
        s.units.length = 'angstrom'

        if len(self.molecule.lattice) in [1,2] and self.molecule.align_lattice():
            log("The lattice supplied for job {} did not follow the convention required by BAND. I rotated the whole system for you. You're welcome".format(self._full_name()), 3)

        for i,atom in enumerate(self.molecule):
            s.atoms['_'+str(i+1)] = atom.str(symbol=self._atom_symbol(atom), space=18, decimal=10)

        if self.molecule.lattice:
            for i,vec in enumerate(self.molecule.lattice):
                s.lattice['_'+str(i+1)] = '{:16.10f} {:16.10f} {:16.10f}'.format(*vec)


    def _remove_mol(self):
        if 'atoms' in self.settings.input:
            del self.settings.input.atoms
        if 'lattice' in self.settings.input:
            del self.settings.input.lattice