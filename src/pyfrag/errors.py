from typing import Optional


class UnsupportedSectionError(Exception):
    """Error raised when the section is not supported when trying to extract sections from the input file."""

    ...


class PyFragSectionInputError(ValueError):
    """An error that occurs when the PyFrag input is invalid."""

    def __init__(self, message, key: Optional[str] = None):
        if key is not None:
            message = f"{key} is not valid. {message}"
        super().__init__(message)
        self.key = key
