import logging
import pprint
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Literal, Union

from pyfrag.config.config import pyfrag_config
from pyfrag.errors import PyFragSectionInputError

logger = logging.getLogger(name="PyFrag Section Reader")


@dataclass
class InputBlocks:
    """Dataclass to store the input blocks of the PyFrag input file."""

    JOBSUB: str = ""
    PYFRAG: str = ""
    ADF: str = ""
    AMS: str = ""
    FRAGMENT_EXTRA: List[str] = field(default_factory=list)
    FRAGMENT_OPEN_EXTRA: List[str] = field(default_factory=list)
    COMPLEX_EXTRA: str = ""

    def get_input_blocks_specific_section(self, section_name: Literal["SCM", "PYFRAG", "JOBSUB"]) -> Dict[str, Union[str, None]]:
        """Returns the content of the specific section. Note: SCM"""
        section_name_str = section_name.upper()

        if section_name_str not in self.__dict__.keys() and section_name_str != "SCM":
            raise PyFragSectionInputError(f"Section {section_name} is not present in the input file.")

        if section_name_str == "SCM":
            ams_keys = pyfrag_config.ams.input_blocks
            return {k: v for k, v in self.__dict__.items() if k in ams_keys}
        else:
            return self.__dict__[section_name.upper()]

    def validate(self):
        if not self.PYFRAG:
            raise PyFragSectionInputError("PyFrag section is missing in the input file.\nPlease make sure to include 'PyFrag' ... ''PyFrag END' in the input file.")

        if not self.ADF and not self.AMS:
            raise PyFragSectionInputError("ADF and AMS sections are missing in the input file.\nPlease make sure to include 'ADF' ... ''ADF END' or 'AMS' ... ''AMS END' in the input file.")

    def __iter__(self):
        for k, v in self.__dict__.items():
            if k == "FRAGMENT_EXTRA":
                for i, fragment in enumerate(v, start=1):
                    yield f"FRAGMENT{i}_EXTRA", fragment
            elif k == "FRAGMENT_OPEN_EXTRA":
                for i, fragment in enumerate(v, start=1):
                    yield f"FRAGMENT{i}_OPEN_EXTRA", fragment
            elif v:
                yield k, v

    def __getitem__(self, key: str):
        return self.__dict__[key]

    def __str__(self) -> str:
        return str(pprint.pformat(asdict(self)))

    def __contains__(self, key: str):
        return key in self.__dict__.keys()
