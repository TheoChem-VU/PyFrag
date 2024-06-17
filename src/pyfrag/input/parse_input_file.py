import difflib
import functools
import re
from pathlib import Path
from typing import Dict, Set, Union

from pyfrag.errors import UnsupportedSectionError

AVAILABLE_SECTIONS = ["JOBSUB", "PYFRAG", "ADF", "AMS", "FRAGMENT1 EXTRA", "FRAGMENT2 EXTRA", "COMPLEX EXTRA"]

_SPELLING_MISTAKES_SECTION_NAMES = {
    "JOBSUB": ["JOB_SUB", "JOB SUB", "JOB_SUBMIT", "JOB SUBMIT"],
    "FRAGMENT1 EXTRA": ["FRAGMENT1_EXTRA"],
    "FRAGMENT2 EXTRA": ["FRAGMENT2_EXTRA"],
    "COMPLEX EXTRA": ["COMPLEX_EXTRA"],
}

ALL_SECTIONS = set(AVAILABLE_SECTIONS + [section for sections in _SPELLING_MISTAKES_SECTION_NAMES.values() for section in sections])


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


def extract_sections(file_content: str) -> Dict[str, str]:
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
        "PYFRAG": [],
        "ADF": [],
        "AMS": [],
        "FRAGMENT1 EXTRA": [],
        "FRAGMENT2 EXTRA": [],
        "COMPLEX EXTRA": [],

    """
    sections_content = {section: "" for section in AVAILABLE_SECTIONS}

    # Extract sections
    for section in AVAILABLE_SECTIONS:
        pattern = f"{section}\n(.*?){section} END"
        match = re.search(pattern, file_content, re.DOTALL | re.IGNORECASE)
        if match:
            sections_content[section] = match.group(1).strip()
        print(f"Extracted {section} section: {sections_content[section]}")
    return sections_content

def main():
    from scm.plams import AMSJob, Molecule, Settings, config, finish, init
    from tcutility.results.ams import get_ams_input
    file_path = Path(__file__).parent.parent.parent.parent / "tests" / "fixtures" / "input_file_examples" / "test.in"
    print(f"Reading file: {file_path}\nExists: {file_path.exists()}")

    with open(file_path, 'r') as f:
        file_content = f.read()

    sections_content = extract_sections(file_content)

    for section, content in sections_content.items():
        if not content:
            continue

        print(f"--- {section} ---")
        print(content)

    # Only get the AMS section as a string

    s = get_ams_input(sections_content["COMPLEX EXTRA"])
    sett = Settings()
    sett.input= s.as_plams_settings()
    print(sett.input)


    correct_frag_sett = Settings()
    correct_frag_sett.input.adf.FragOccupations
    correct_frag_sett.input.adf.FragOccupations._1 = "CX"
    correct_frag_sett.input.adf.FragOccupations._2 = "A1 10 // 9"
    correct_frag_sett.input.adf.FragOccupations._3 = "A2 0 // 0"
    correct_frag_sett.input.adf.FragOccupations._4 = "B1 3 // 3"
    correct_frag_sett.input.adf.FragOccupations._5 = "B2 5 // 4"
    correct_frag_sett.input.adf.FragOccupations._6 = "subend"
    correct_frag_sett.input.adf.FragOccupations._7 = "NH2"
    correct_frag_sett.input.adf.FragOccupations._8 = "A1 3 // 4"
    correct_frag_sett.input.adf.FragOccupations._9 = "A2 1 // 1"
    correct_frag_sett.input.adf.FragOccupations._10 = "B1 1 // 1"
    correct_frag_sett.input.adf.FragOccupations._11 = "B2 3 // 4"
    correct_frag_sett.input.adf.FragOccupations._12 = "subend"


    correct_frag_sett.input.adf.IrrepOccupations.A1 = "18 1"
    correct_frag_sett.input.adf.IrrepOccupations.A2 = "0"
    correct_frag_sett.input.adf.IrrepOccupations.B1 = "6"
    correct_frag_sett.input.adf.IrrepOccupations.B2 = "8 1"
    print(correct_frag_sett)


    non_standard_blocks = ["fragoccupations", "irrepoccupations"]
    # Iterate over the fragoccupations list
    for key in sett.input.adf:
        if key not in non_standard_blocks:
            continue

        for i, occupation in enumerate(sett.input.adf.pop(key), start=1):
            # Add each occupation as a key-value pair to the FragOccupations dictionary
            sett.input.adf[key]['_' + str(i)] = occupation


    init()
    config.erase_workdir = True
    config.preview = True
    mol = Molecule.from_elements("OO")
    job = AMSJob(settings = sett, name = "test", molecule = mol)

    print(sett.input.adf)
    print(correct_frag_sett.input.adf)
    print(job.get_input())
    finish()

if __name__ == "__main__":
    main()
