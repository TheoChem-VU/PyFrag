import pathlib as pl
import shutil
from abc import ABC, abstractmethod
from typing import List, Union

from scm.plams import MultiJob, config, init, load_all

from pyfrag.enums import JobStatus
from pyfrag.jobs.fragment import Complex, Fragment
from pyfrag.jobs.fragmentjob import PyFragFragmentJob


def _check_for_restart_folder(job_dir: pl.Path) -> pl.Path:
    res_dir = job_dir.with_suffix(".res")
    if job_dir.exists():
        if res_dir.exists():
            shutil.rmtree(res_dir)
        job_dir.rename(res_dir)
    return res_dir


def _remove_restart_folder():
    job_dir = pl.Path(config.default_jobmanager.workdir)
    res_dir = job_dir.with_suffix(".res")
    if res_dir.exists():
        shutil.rmtree(res_dir)


def _setup_plams_environment(job_name: str, job_dir: pl.Path):
    if job_dir.exists():
        res_dir = _check_for_restart_folder(job_dir / job_name)

    init(path=str(job_dir), folder=job_name, use_existing_folder=False)

    load_all(res_dir) if res_dir.exists() else None

    config.job.keep = ["-", "$JN.err", "CreateAtoms.out", "t12.rel"]  # type: ignore


class PyFragBaseJob(ABC):
    def __init__(self, job_dir: pl.Path, complex: Complex, fragments: List[Fragment]):
        self._multi_jobs: Union[list[MultiJob], None] = None
        self.complex = complex
        self.fragments = fragments
        self.job_dir = job_dir
        self.job_status: JobStatus = JobStatus.INITIALIZED

    def __enter__(self, *args, **kwargs):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        _remove_restart_folder()
        self.finish()

    @property
    def len_trajectory(self):
        return len(self.complex.trajectory)

    @abstractmethod
    def prepare_job(self):
        pass

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def finish(self):
        pass


class RestrictedPyFragJob(PyFragBaseJob):
    def prepare_job(self):
        # Create all multijobs
        self._multi_jobs = []

        for i in range(self.len_trajectory):
            frag_mols = [frag.trajectory[i] for frag in self.fragments]
            frag_settings = [frag.plams_settings for frag in self.fragments]

            pyfrag_job = PyFragFragmentJob(
                name=f"{i+1}-{self.len_trajectory}",
                full_settings=self.complex.plams_settings,
                frag_mols=frag_mols,
                frag_settings=frag_settings,
                fragment_names=[frag.name for frag in self.fragments],
                complex_name=self.complex.name,
            )

            self._multi_jobs.append(pyfrag_job)

        self.job_status = JobStatus.PREPARED

    def run(self):
        # Initialize the PLAMS environment
        _setup_plams_environment(job_name="pyfrag_job", job_dir=self.job_dir)

        if not self.job_status == JobStatus.PREPARED:
            self.prepare_job()

        for job in self._multi_jobs:  # type: ignore # this is always an instantiated list of multijobs!
            job.run(jobmanager=config.default_jobmanager)

        self.job_status = JobStatus.RUNNING

    def finish(self):
        pass
