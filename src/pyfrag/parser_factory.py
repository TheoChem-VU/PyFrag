from enum import StrEnum
from typing import Callable, Dict

import pyfrag.executables.adf.adf_parser as adf_parser


class ExecutableType(StrEnum):
    ADF = "adf"


ParserFunc = Callable[[str], Dict[str, str]]


def get_parser_from_executable(executable: str) -> ParserFunc:
    """Get the parser for a specific PyFrag executable."""

    if executable.lower() == ExecutableType.ADF:
        return adf_parser.extract_sections
    else:
        raise ValueError(f"Unsupported executable: {executable}. Supported: {list(ExecutableType)}.")
