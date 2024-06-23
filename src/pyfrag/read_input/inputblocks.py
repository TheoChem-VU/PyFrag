import logging
from dataclasses import asdict, dataclass
from pprint import pprint
from typing import Union

from pyfrag.errors import PyFragSectionInputError

logger = logging.getLogger(name="PyFrag Section Reader")


@dataclass
class InputBlocks:
    """Dataclass to store the input blocks of the PyFrag input file."""

    JOBSUB: Union[str, None] = None
    PYFRAG: Union[str, None] = None
    ADF: Union[str, None] = None
    AMS: Union[str, None] = None
    FRAGMENT1_EXTRA: Union[str, None] = None
    FRAGMENT2_EXTRA: Union[str, None] = None
    COMPLEX_EXTRA: Union[str, None] = None

    def validate(self):
        if self.PYFRAG is None:
            raise PyFragSectionInputError("PyFrag section is missing in the input file.\nPlease make sure to include 'PyFrag' ... ''PyFrag END' in the input file.")

        if self.ADF is None and self.AMS is None:
            raise PyFragSectionInputError("ADF and AMS sections are missing in the input file.\nPlease make sure to include 'ADF' ... ''ADF END' or 'AMS' ... ''AMS END' in the input file.")

    def __iter__(self):
        return iter((k, v) for k, v in self.__dict__.items() if v is not None)

    def __getitem__(self, key: str):
        return self.__dict__[key]

    def __str__(self) -> str:
        return str(pprint(asdict(self)))
