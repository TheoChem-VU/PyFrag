from enum import StrEnum
from typing import Callable, Dict

import pyfrag.executables.adf.adf_parser as adf_parser


class ExecutableType(StrEnum):
    ADF = "adf"
    ORCA = "orca"
    GAUSSIAN = "gaussian"
    TURBOMOLE = "turbomole"


ParserFunc = Callable[[str], Dict[str, str]]


def get_parser() -> ParserFunc:
    """
        Get the parser for a specific PyFrag executable.
        Currently, all the executables use the same parser. It may be extended in the future.

    The `adf_parser.extract_sections` function is used by the `main.py` script to extract the `jobsub` section that specifies the job submission parameters (to a job scheduler such as SLURM or PBS).
    """
    return adf_parser.extract_sections
