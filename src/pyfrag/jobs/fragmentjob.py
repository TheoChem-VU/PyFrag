from typing import List, Optional

from scm.plams import AMSJob, Molecule, MultiJob, Settings


class PyFragFragmentJob(MultiJob):
    def __init__(
        self,
        full_settings: Settings,
        frag_mols: List[Molecule],
        frag_settings: Optional[List[Settings]] = None,
        fragment_names: Optional[List[str]] = None,
        complex_name: str = "full",
        **kwargs,
    ):
        MultiJob.__init__(self, **kwargs)
        self.full_settings = full_settings or Settings()
        self.frag_mols = [frag.copy() for frag in frag_mols]
        self.frag_settings = frag_settings or [Settings() for _ in range(len(frag_mols))]
        self.fragment_names = [f"f{i+1}" for i in range(len(frag_mols))] if fragment_names is None else fragment_names
        self.complex_name: str = complex_name
        self.frags: List[AMSJob] = list()
        self.children: List[AMSJob] = list()

    def prerun(self):
        for frag_mol, frag_settings, frag_name in zip(self.frag_mols, self.frag_settings, self.fragment_names):
            # Please note that the settings of the fragments ARE LEADING over the full settings (soft update)
            self.frags.append(AMSJob(name=frag_name, molecule=frag_mol, settings=frag_settings.soft_update(self.settings)))

        for frag_name, frag_mol in zip(self.fragment_names, self.frag_mols):
            for at in frag_mol:
                at.properties.suffix = f"adf.f={frag_name}"

        complex_mol = sum(self.frag_mols, Molecule())
        self.full = AMSJob(name=self.complex_name, molecule=complex_mol, settings=self.full_settings.soft_update(self.settings))

        for frag_name, frag_job in zip(self.fragment_names, self.frags):
            self.full.settings.input.adf.fragments[frag_name] = (frag_job, "adf")

        self.children = self.frags + [self.full]
