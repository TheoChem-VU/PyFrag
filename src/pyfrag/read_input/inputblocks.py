import logging
import pprint
from dataclasses import asdict, dataclass, field
from typing import Dict, Literal, Union

from pyfrag.errors import PyFragSectionInputError

logger = logging.getLogger(name="PyFrag Section Reader")


@dataclass
class InputBlocks:
    """Dataclass that stores the input blocks of the PyFrag input file."""

    JOBSUB: str = ""
    PYFRAG: str = ""
    ADF: str = ""
    AMS: str = ""
    FRAGMENT_EXTRA: Dict[int, str] = field(default_factory=dict)
    FRAGMENT_OPEN_EXTRA: Dict[int, str] = field(default_factory=dict)
    COMPLEX_EXTRA: str = ""

    def get_input_blocks_specific_section(self, section_name: Literal["SCM", "PYFRAG", "JOBSUB"]) -> Dict[str, Union[str, None]]:
        """Returns the content of the specific section. Note: SCM includes both the PyFrag-defined AMS and ADF sections (AMS is the >2019 parsing of SCM and ADF is the parsing format before 2019)."""
        section_name_str = section_name.upper()

        if section_name_str not in self.__dict__.keys() and section_name_str != "SCM":
            raise PyFragSectionInputError(f"Section {section_name} is not present in the input file.")

        if section_name_str == "SCM":
            ams_keys = ["AMS", "ADF"]
            return {k: v for k, v in self.__dict__.items() if k in ams_keys}
        else:
            return self.__dict__[section_name.upper()]

    def validate(self):
        if not self.PYFRAG:
            raise PyFragSectionInputError("PyFrag section is missing in the input file.\nPlease make sure to include 'PyFrag' ... ''PyFrag END' in the input file.")

        if not self.ADF and not self.AMS:
            raise PyFragSectionInputError("ADF and AMS sections are missing in the input file.\nPlease make sure to include 'ADF' ... ''ADF END' or 'AMS' ... ''AMS END' in the input file.")

        if self.ADF and self.AMS:
            raise PyFragSectionInputError("Both ADF and AMS sections are present in the input file.\nPlease only include one of the two sections, preferably the AMS section.")

    def __iter__(self):
        for k, v in self.__dict__.items():
            if k == "FRAGMENT_EXTRA":
                for i, fragment in v.items():
                    yield f"FRAGMENT{i}_EXTRA", fragment
            elif k == "FRAGMENT_OPEN_EXTRA":
                for i, fragment in v.items():
                    yield f"FRAGMENT{i}_OPEN_EXTRA", fragment
            elif v:
                yield k, v

    def __getitem__(self, key: str):
        return self.__dict__[key]

    def __str__(self) -> str:
        return str(pprint.pformat(asdict(self)))

    def __contains__(self, key: str):
        return key in self.__dict__.keys()
