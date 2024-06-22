class UnsupportedSectionError(Exception):
    """Error raised when the section is not supported when trying to extract sections from the input file."""

    ...


class PyFragSectionInputError(Exception):
    """Error raised when the PyFrag section is not found in the input file."""

    ...
