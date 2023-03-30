from .scmjob import SCMJob, SCMResults

__all__ = ['UFFJob', 'UFFResults']



class UFFResults(SCMResults):
    _kfext = '.rkf'
    _rename_map = {'uff.rkf':'$JN'+_kfext}



class UFFJob(SCMJob):
    _result_type = UFFResults
    _command = 'uff'
    _top = ['title', 'units']
    _subblock_end = 'end'

    def _serialize_mol(self):
        s = self.settings.input

        printtypes = all(map(lambda at: ('uff' in at.properties and 'type' in at.properties.uff), self.molecule))

        for i,atom in enumerate(self.molecule):
            suffix = ''
            suffix_dict = {}
            if printtypes:
                suffix = '  {type}  {charge}'
                suffix_dict = atom.properties.uff.copy()
                if 'charge' not in suffix_dict:
                    suffix_dict.charge = 0.0
            s.system.atoms['_'+str(i+1)] = atom.str(suffix=suffix, suffix_dict=suffix_dict, space=18, decimal=10)

        if 'ignore_bonds' not in self.settings:
            self.molecule.set_atoms_id()
            for i,bond in enumerate(self.molecule.bonds):
                s.system.bonds['_'+str(i+1)] = '{:6d}  {:6d}  {:6.2f}'.format(bond.atom1.id, bond.atom2.id, bond.order)
            self.molecule.unset_atoms_id()

        if self.molecule.lattice:
            for i,vec in enumerate(self.molecule.lattice):
                s.system.lattice['_'+str(i+1)] = '{:16.10f} {:16.10f} {:16.10f}'.format(*vec)

    def _remove_mol(self):
        s = self.settings.input

        if 'system' in s:
            for k in ['atoms', 'bonds', 'lattice']:
                if k in s.system:
                    del s.system[k]
