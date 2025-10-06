from typing import List, Optional

# =============================================================================
# Base PyFrag error
# =============================================================================


class PyFragError(Exception):
    """Base class for PyFrag errors."""


# =============================================================================
# Unsupported executable error
# =============================================================================


class ExecutablePathNotFoundError(PyFragError):
    """Error raised when the executable path is not found."""


class ExecutableNotSupportedError(PyFragError):
    """Error raised when an unsupported executable is requested."""

    def __init__(self, executable: str, valid_executables: List[str]):
        self.executable = executable
        self.valid_executables = valid_executables
        super().__init__(f"Executable '{executable}' is not supported. Valid options are: {', '.join(valid_executables)}")


# =============================================================================
# AMS environment error
# =============================================================================


class AMSNotFoundError(PyFragError):
    """Error raised when the AMS environment is not found."""


# =============================================================================
# Input reading and parsing errors
# =============================================================================


class PyFragInputFileNotFoundError(PyFragError):
    """Error raised when the input file is not found."""


class UnsupportedSectionError(PyFragError):
    """Error raised when the section is not supported when trying to extract sections from the input file."""


class PyFragSectionInputError(PyFragError):
    """An error that occurs when the PyFrag input is invalid."""

    def __init__(self, message, key: Optional[str] = None):
        if key is not None:
            message = f"{key} is not valid. {message}"
        super().__init__(message)
        self.key = key


class AMSSectionInputError(PyFragError):
    """An error that occurs when the PyFrag input is invalid."""

    def __init__(self, message, key: Optional[str] = None):
        if key is not None:
            message = f"{key} is not valid. {message}"
        super().__init__(message)
        self.key = key


# =============================================================================
# Coordinate file errors
# =============================================================================
class PyFragCoordFileError(PyFragError):
    """An error that occurs when the coordinate file is invalid."""


class FragmentIndicesError(PyFragError):
    """An error that occurs when the fragment indices are invalid."""


# =============================================================================
# AMS input error
# =============================================================================


class PyFragAMSSectionReadError(PyFragError):
    """An error that occurs when the AMS section is invalid."""


# =============================================================================
# Fragment optimization error
# =============================================================================


class FragmentOptimizationError(PyFragError):
    """An error that occurs when the fragment optimization fails."""


# =============================================================================
# PyFrag job driver error
# =============================================================================


class PyFragSortComplexMoleculeError(PyFragError):
    """An error that occurs when the PyFrag driver fails when sorting the complex molecule to the original fragment input order"""
