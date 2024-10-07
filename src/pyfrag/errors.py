from typing import Optional

# =============================================================================
# Base PyFrag error
# =============================================================================


class PyFragError(Exception):
    """Base class for PyFrag errors."""


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


# =============================================================================
# AMS input error
# =============================================================================


class PyFragAMSSectionReadError(PyFragError):
    """An error that occurs when the AMS section is invalid."""
