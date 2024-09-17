import pathlib as pl
from abc import ABC, abstractmethod
from typing import List, Union

from scm.plams import MultiJob, config, init


def _setup_plams_environment(job_name: str, job_dir: pl.Path):
    init(path=str(job_dir), folder=job_name, use_existing_folder=True)

    config.job.keep = ["-", "$JN.err", "CreateAtoms.out", "t12.rel"]  # type: ignore


class PyFragBaseJob(ABC):
    def __init__(self, trajectories: List[List[int]], **kwargs):
        self.multi_jobs: Union[list[MultiJob], None] = None
        self.molecule_trajectories: List[List[int]] = trajectories

    def __enter__(self, job_name: str = "pyfrag_job", job_dir: pl.Path = pl.Path.cwd()):
        _setup_plams_environment(job_name, job_dir)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finish()

    @abstractmethod
    def prepare_job(self):
        pass

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def finish(self):
        pass
