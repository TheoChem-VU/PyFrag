from typing import Optional


class PyFragError(Exception):
    """Base class for PyFrag errors."""


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


class PyFragCoordFileError(PyFragError):
    """An error that occurs when the coordinate file is invalid."""


class PyFragAMSSectionReadError(PyFragError):
    """An error that occurs when the AMS section is invalid."""
