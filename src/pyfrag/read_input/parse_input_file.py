import logging
import re
from pathlib import Path
from typing import Union

from pyfrag.errors import PyFragInputFileNotFoundError
from pyfrag.read_input.inputblocks import InputBlocks

logger = logging.getLogger(name="PyFrag Section Reader")


# =============================================================
# File Handling Functions =====================================
# =============================================================


def read_input_file_content(file_path: Union[str, Path]) -> str:
    """Reads the content of a given file. May raises errors when the file does not exist."""
    file = Path(file_path)

    if not file.exists():
        raise PyFragInputFileNotFoundError(f"Input file {file_path} does not exist. Please make sure the file exists.")

    with open(file, "r") as f:
        return f.read()


# =============================================================
# Section Reading Functions ===================================
# =============================================================


def correct_spelling_mistakes(file_content: str) -> str:
    """
    Correct spelling mistakes in the input file content. Currently, it tests for the following cases:

        - FRAGMENT1_EXTRA -> FRAGMENT1 EXTRA
        - FRAGMENT1_OPEN_EXTRA -> FRAGMENT1 OPEN EXTRA
        - COMPLEX_EXTRA -> COMPLEX EXTRA
        - JOB SUB -> JOBSUB
        - JOB_SUB -> JOBSUB

    """
    # Correct section names like FRAGMENT1_EXTRA to FRAGMENT1 EXTRA and FRAGMENT1_OPEN_EXTRA to FRAGMENT1 OPEN EXTRA.
    pattern_fragment = re.compile(r"(FRAGMENT\d+)_(OPEN_)?EXTRA", re.IGNORECASE)
    corrected_content = pattern_fragment.sub(r"\1 \2EXTRA", file_content)

    # Correct section names like JOB SUB and JOB_SUB to JOBSUB.
    pattern_jobsub = re.compile(r"JOB[\s_]+SUB", re.IGNORECASE)
    corrected_content = pattern_jobsub.sub("JOBSUB", corrected_content)

    return corrected_content


def remove_comments(section_content: str, section_name: str) -> str:
    """Removes comments from the input file content which are lines starting with # or !."""

    # JOBSUB contains sbatch options such as #SBATCH -N 1, which should not be removed. Thus, return the content as is.
    if section_name == "JOBSUB":
        return section_content

    lines = section_content.split("\n")
    section_content_without_comments = [line for line in lines if not line.strip().startswith(("#", "!"))]
    return "\n".join(section_content_without_comments)


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
        "FRAGMENT EXTRA": [" "],
        "COMPLEX EXTRA": "",

    Args:
        file_content (str): The content of the input file.

    Returns:
        InputBlocks: A dictionary with the section name as key and the content as value.


    """
    sections_content = InputBlocks()

    # Preprocess the file content
    preprocessed_content = correct_spelling_mistakes(file_content)

    # This pattern matches the section name (case-insensitive) followed by the content of the section. There is no limit to how many "FRAGMENT{x} EXTRA" sections there can be.
    # Note that this pattern also graps sections that are commented out. A fix is hard because the "JOBSUB" section contains sbatch lines that start also with a comment
    pattern = re.compile(r"(JOBSUB|PYFRAG|ADF|AMS|COMPLEX EXTRA|FRAGMENT\d+ EXTRA|FRAGMENT\d+ OPEN EXTRA)\n(.*?)(\1 END)", re.DOTALL | re.IGNORECASE)

    # Find all matches
    matches = pattern.findall(preprocessed_content)

    for match in matches:
        section_name = match[0].upper()
        section_content = match[1].strip()

        section_content = remove_comments(section_content, section_name)

        if section_name == "JOBSUB":
            sections_content.JOBSUB = section_content
        elif section_name == "PYFRAG":
            sections_content.PYFRAG = section_content
        elif section_name == "ADF":
            sections_content.ADF = section_content
        elif section_name == "AMS":
            sections_content.AMS = section_content
        elif section_name == "COMPLEX EXTRA":
            sections_content.COMPLEX_EXTRA = section_content
        elif re.match(r"FRAGMENT\d+ EXTRA", section_name, re.IGNORECASE):
            sections_content.FRAGMENT_EXTRA.append(section_content)
        elif re.match(r"FRAGMENT\d+ OPEN EXTRA", section_name, re.IGNORECASE):
            sections_content.FRAGMENT_OPEN_EXTRA.append(section_content)

    return sections_content
