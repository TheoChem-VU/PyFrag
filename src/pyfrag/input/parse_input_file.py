import difflib
import functools
import logging
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from pprint import pprint
from typing import Set, Union

from pyfrag.errors import UnsupportedSectionError

logger = logging.getLogger()

AVAILABLE_SECTIONS = ["JOBSUB", "PYFRAG", "ADF", "AMS", "FRAGMENT1 EXTRA", "FRAGMENT2 EXTRA", "COMPLEX EXTRA"]

_SPELLING_MISTAKES_SECTION_NAMES = {
    "JOBSUB": ["JOB_SUB", "JOB SUB", "JOB_SUBMIT", "JOB SUBMIT"],
    "FRAGMENT1 EXTRA": ["FRAGMENT1_EXTRA"],
    "FRAGMENT2 EXTRA": ["FRAGMENT2_EXTRA"],
    "COMPLEX EXTRA": ["COMPLEX_EXTRA"],
}

ALL_SECTIONS = set(AVAILABLE_SECTIONS + [section for sections in _SPELLING_MISTAKES_SECTION_NAMES.values() for section in sections])


@dataclass
class InputBlocks:
    JOBSUB: Union[str, None] = None
    PYFRAG: Union[str, None] = None
    ADF: Union[str, None] = None
    AMS: Union[str, None] = None
    FRAGMENT1_EXTRA: Union[str, None] = None
    FRAGMENT2_EXTRA: Union[str, None] = None
    COMPLEX_EXTRA: Union[str, None] = None

    def __iter__(self):
        return iter((k, v) for k, v in self.__dict__.items() if v is not None)

    def __getitem__(self, key: str):
        return self.__dict__[key]

    def __str__(self) -> str:
        return str(pprint(asdict(self)))


@functools.lru_cache(1)
def _get_all_sections() -> Set[str]:
    set_sections = set(AVAILABLE_SECTIONS)
    for sections in _SPELLING_MISTAKES_SECTION_NAMES.values():
        set_sections.update(sections)
    return set_sections


def get_one_section(line_upper: str) -> Union[str, None]:
    """Get the section name from the line. If the section is not supported, raise an error and suggest the closest supported section."""
    all_sections = _get_all_sections()
    close_matches = difflib.get_close_matches(line_upper, all_sections, n=1, cutoff=0.8)
    if close_matches:
        current_section = close_matches[0]
        if current_section not in AVAILABLE_SECTIONS:
            raise UnsupportedSectionError(f"Unsupported section: {current_section}. Please use one of {' / '.join(difflib.get_close_matches(current_section, AVAILABLE_SECTIONS))} instead.")
        return current_section


def extract_section_blocks_from_file_content(file_content: str) -> InputBlocks:
    r"""
    Extract the content of each section from the input file.
    Return a dictionary with the section name as key and the content as value.

    A section is given by the following format:
    SECTION_NAME
        content
        content
        ...
    SECTION_NAME END

    Example:
    JOBSUB
        #!/bin/bash
        #SBATCH -N 1
    JOBSUB END

    The function will return:
    {
        "JOBSUB": ["#!/bin/bash\n", "#SBATCH -N 1\n"],
        "PYFRAG": "",
        "ADF": "",
        "AMS": "",
        "FRAGMENT1 EXTRA": "",
        "FRAGMENT2 EXTRA": "",
        "COMPLEX EXTRA": "",

    """
    sections_content: InputBlocks = InputBlocks()
    # Extract sections
    for section in AVAILABLE_SECTIONS:
        pattern = f"{section}\n(.*?){section} END"
        match = re.search(pattern, file_content, re.DOTALL | re.IGNORECASE)
        if match:
            sections_content.__setattr__(section, match.group(1).strip())
    return sections_content


def main():
    from scm.plams import AMSJob, Molecule, Settings, config, finish, init
    from tcutility.results.ams import get_ams_input

    file_path = Path(__file__).parent.parent.parent.parent / "tests" / "fixtures" / "input_file_examples" / "test.in"
    print(f"Reading file: {file_path}\nExists: {file_path.exists()}")

    with open(file_path, "r") as f:
        file_content = f.read()

    sections_content = extract_section_blocks_from_file_content(file_content)

    s = get_ams_input(sections_content.AMS) if sections_content.AMS else Settings()

    sett = Settings()
    sett.input = s.as_plams_settings()
    print(sett.input)

    correct_frag_sett = Settings()
    correct_frag_sett.input.adf.FragOccupations = {
        "_1": "CX",
        "_2": "A1 10 // 9",
        "_3": "A2 0 // 0",
        "_4": "B1 3 // 3",
        "_5": "B2 5 // 4",
        "_6": "subend",
        "_7": "NH2",
        "_8": "A1 3 // 4",
        "_9": "A2 1 // 1",
        "_10": "B1 1 // 1",
        "_11": "B2 3 // 4",
        "_12": "subend",
    }

    correct_frag_sett.input.adf.IrrepOccupations = {"A1": "18 1", "A2": "0", "B1": "6", "B2": "8 1"}
    print(correct_frag_sett)

    non_standard_blocks = ["fragoccupations", "irrepoccupations"]
    # Iterate over the fragoccupations list
    for key in sett.input.adf:
        if key not in non_standard_blocks:
            continue

        for i, occupation in enumerate(sett.input.adf.pop(key), start=1):
            # Add each occupation as a key-value pair to the FragOccupations dictionary
            sett.input.adf[key]["_" + str(i)] = occupation

    init()
    config.erase_workdir = True
    config.preview = True
    mol = Molecule.from_elements("OO")
    job = AMSJob(settings=sett, name="test", molecule=mol)

    print(sett.input.adf)
    print(correct_frag_sett.input.adf)
    print(job.get_input())
    finish()


if __name__ == "__main__":
    main()
