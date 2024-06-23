import functools
import logging
import re
from pathlib import Path
from typing import Set, Union

from pyfrag.input.inputblocks import InputBlocks

logger = logging.getLogger(name="PyFrag Section Reader")

AVAILABLE_SECTIONS = ["JOBSUB", "PYFRAG", "ADF", "AMS", "FRAGMENT1 EXTRA", "FRAGMENT2 EXTRA", "COMPLEX EXTRA"]

_SPELLING_MISTAKES_SECTION_NAMES = {
    "JOB SUB": "JOBSUB",
    "JOB_SUB": "JOBSUB",
    "FRAGMENT1_EXTRA": "FRAGMENT1 EXTRA",
    "FRAGMENT2_EXTRA": "FRAGMENT2 EXTRA",
    "COMPLEX_EXTRA": "COMPLEX EXTRA",
}

# =============================================================
# File Handling Functions =====================================
# =============================================================


def read_input_file_content(file_path: Union[str, Path]) -> str:
    """Reads the content of a given file. May raises errors when the file does not exist."""
    file = Path(file_path)

    if not file.exists():
        raise FileNotFoundError(f"File {file_path} does not exist. Please make sure the file exists.")

    with open(file, "r") as f:
        return f.read()


# =============================================================
# Section Reading Functions ===================================
# =============================================================


@functools.lru_cache(1)
def _get_all_sections() -> Set[str]:
    set_sections = set(AVAILABLE_SECTIONS)
    for sections in _SPELLING_MISTAKES_SECTION_NAMES.keys():
        set_sections.add(sections)
    return set_sections


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
    for section in _get_all_sections():
        pattern = f"{section}\n(.*?){section} END"
        match = re.search(pattern, file_content, re.DOTALL | re.IGNORECASE)
        if match and section in AVAILABLE_SECTIONS:
            sections_content.__setattr__(section, match.group(1).strip())
        elif match and section in _SPELLING_MISTAKES_SECTION_NAMES:
            logger.warning(f"Spelling mistake found in section name: {section}. Correcting it to {_SPELLING_MISTAKES_SECTION_NAMES[section]}.")
            sections_content.__setattr__(_SPELLING_MISTAKES_SECTION_NAMES[section], match.group(1).strip())
    return sections_content
