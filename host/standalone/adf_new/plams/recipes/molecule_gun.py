from ..interfaces.adfsuite.reaxff import ReaxFFJob

__all__ = ['MoleculeGunJob']

class MoleculeGunJob(ReaxFFJob):

    def __init__(self, bullet, **kwargs):
        ReaxFFJob.__init__(self, **kwargs)
        self.bullet = bullet

    def _get_ready(self):
        ReaxFFJob._get_ready(self)
        self._write_geofile(molecule=self.bullet, filename='addmol.bgf', settings=self.settings.input.molecule_gun, description=self.bullet.get_formula(), lattice=False)
