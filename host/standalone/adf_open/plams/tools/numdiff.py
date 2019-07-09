from __future__ import unicode_literals

from collections import OrderedDict

from ..core.results import Results
from ..core.basejob import MultiJob
from ..core.settings import Settings
from ..core.basemol import Atom
from ..interfaces.adfsuite import ADFJob, BANDJob, DFTBJob


__all__ = ['ADFNumGradJob' ,'BANDNumGradJob', 'DFTBNumGradJob']

class NumGradResults(Results):

    def get_gradient(self, atom, coord, func=None):
    #func should be a function that takes results of single point job and returns a number. If None, get_energy is used.

        s = self.job.settings.numgrad
        extract = func or s.get_energy
        energies = [extract(self.job.children[(atom,coord,i)].results) for i in s.steps]
        coeffs, denom = self.job._coeffs[s.npoints]
        return sum([c*e for c,e in zip(coeffs, energies)])/(denom * s.step)


class NumGradJob(MultiJob):
    _result_type = NumGradResults
    _coeffs = {
        3: ([-1,1], 2),
        5: ([1,-8,8,-1], 12),
        7: ([-1,9,-75,75,-9,1], 60),
        9: ([3,-32,168,-672,672,-168,32,-3], 840),
    }

    def __init__(self, molecule, atoms=None, npoints=3, step=0.01, unit='angstrom', **kwargs):
        #atoms should be a list. Each element is the id of atom or tuple like (2,'xz') if gradients along not all direction are to be calculated. Instead of atom's id the pointer to actual Atom object can be passed.

        MultiJob.__init__(self, children=OrderedDict(), **kwargs)
        s = self.settings.numgrad
        s.molecule = molecule
        s.atoms = atoms
        s.npoints = npoints
        s.step = step
        s.unit = unit

        #this is an abstract class, the following two guys have to be set in subclasses
        s.jobtype = None   #who is going to do the actual work and calculate single points
        s.get_energy = None   #should be a function that takes results of single point job and returns a number

    def prerun(self):
        s = self.settings.numgrad

        s.steps = list(range(-(s.npoints//2), 0)) + list(range(1, s.npoints//2+1))

        atomlist = s.atoms or list(range(1, len(s.molecule)+1))
        for atom in atomlist:
            if isinstance(atom, tuple):
                at, co = atom
            else:
                at, co = atom, 'xyz'
            if isinstance(at, Atom):
                at = s.molecule.atoms.index(at)+1
            #now at is id of atom to move and co contains directions

            for axis in co:
                for i in s.steps:
                    v = (s.step*i if axis=='x' else 0, s.step*i if axis=='y' else 0, s.step*i if axis=='z' else 0)
                    newmol = s.molecule.copy()
                    newmol.atoms[at-1].translate(v, s.unit)
                    newname = self.name + str(at) + axis + str(i)
                    self.children[(at,axis,i)] = s.jobtype(name=newname, molecule=newmol, settings=self.settings.child)

#===================================================================================================


def _bond_energy(results):
    return results.readkf('Energy', 'Bond Energy')

class ADFNumGradJob(NumGradJob):
    def __init__(self, **kwargs):
        NumGradJob.__init__(self, **kwargs)
        self.settings.numgrad.jobtype = ADFJob
        self.settings.numgrad.get_energy = _bond_energy

#===================================================================================================

def _BAND_totalenergy(results):
    return results.readkf('Bond energies', 'final bond energy')

class BANDNumGradJob(NumGradJob):
    def __init__(self, **kwargs):
        NumGradJob.__init__(self, **kwargs)
        self.settings.numgrad.jobtype = BANDJob
        self.settings.numgrad.get_energy = _BAND_totalenergy

#===================================================================================================

def _DFTB_totalenergy(results):
    return float(results.grep_output('Total Energy (hartree)')[0].split()[-1])

class DFTBNumGradJob(NumGradJob):
    def __init__(self, **kwargs):
        NumGradJob.__init__(self, **kwargs)
        self.settings.numgrad.jobtype = DFTBJob
        self.settings.numgrad.get_energy = _DFTB_totalenergy

