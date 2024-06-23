"""
Module responsible for parsing the PyFrag input file by splitting the file into sections and reading the sections.
Then, it writes an PyFrag run script (for sbatch purposes) and gives the file path to the main program.
"""

import argparse
import logging
import os
from contextlib import contextmanager
from pathlib import Path

import pyfrag.main_program as pyfrag_driver
from pyfrag.read_input.parse_input_file import InputBlocks, extract_section_blocks_from_file_content, read_input_file_content


def _write_pyfrag_input_file(file_path: Path, input_blocks: InputBlocks) -> None:
    """
    Writes the PyFrag input file to the file path given. It writes the PyFrag section first and then the ADF and AMS sections.
    """
    with open(file_path, "w") as f:
        [f.write(f"\n{section}\n{content}\n{section} END\n") for section, content in input_blocks if section != "JOBSUB"]


def _write_pyfrag_run_script(file_path: Path, input_file_path: Path, input_blocks: InputBlocks) -> None:
    """
    Writes the PyFrag run script to the file path given. First it decides whether to use sbatch or not which depends on the job submission section.
    Then, it writes the sbatch script to the file path.
    """
    jobsub_section = input_blocks.JOBSUB

    with open(file_path, "w") as f:
        f.write(f"{jobsub_section}\n") if jobsub_section is not None else f.write("#!/bin/bash\n")
        f.write(f"\npython {pyfrag_driver.__file__} -f {str(input_file_path)}\n")


@contextmanager
def _run_pyfrag_program(run_file: Path, input_file: Path, use_sbatch: bool):
    """Runs the PyFrag program with the given run file. If use_sbatch is True, it will run the program with sbatch."""
    submit_command = f"sbatch {str(run_file)}" if use_sbatch else f"bash {str(run_file)}"

    try:  # Run the PyFrag program
        os.system(submit_command)
        yield
    finally:
        run_file.unlink(missing_ok=True)


def main():
    """
    This function is responsible for parsing the user input file, making an runscript and submitting/running the PyFrag program.
    1. Read (and validate) the input file content.
    2. Extract the sections from the input file content.
    3. Validate the extracted sections.
    4. Write the PyFrag run script and PyFrag-specific input file.
    5. Run the PyFrag program using a context manager that always removes the run file after the program is done.

    """
    parser = argparse.ArgumentParser(description="Parse PyFrag input file.")
    parser.add_argument("-f", "--file", help="Path to the run file.", required=True)
    parser.add_argument("-v", "--verbose", help="Set the logger to debug mode.", action="store_true")
    args = parser.parse_args()

    file_path = Path(args.file).resolve()
    log_level = logging.getLevelName("DEBUG") if args.verbose else logging.INFO

    logging.basicConfig(level=log_level, format="[%(asctime)s][%(levelname)s][%(module)s]: %(message)s", datefmt="%I:%M:%S")
    logger = logging.getLogger("pyfrag.parser")

    logger.debug(f"Reading file: {file_path.name} (Exists: {file_path.exists()})")
    file_content = read_input_file_content(file_path=file_path)

    section_blocks = extract_section_blocks_from_file_content(file_content)
    section_blocks.validate()
    logger.debug(f"Extracted sections: {', '.join([section for section, _ in section_blocks])}")

    run_file_path = file_path.parent / f"{file_path.stem}.run"
    use_sbatch = section_blocks.JOBSUB is not None

    input_file_path = file_path.parent / f"pyfrag_{file_path.stem}.in"
    _write_pyfrag_input_file(input_file_path, section_blocks)
    _write_pyfrag_run_script(run_file_path, input_file_path, section_blocks)

    with _run_pyfrag_program(run_file_path, use_sbatch=use_sbatch, input_file=input_file_path):
        logger.info(f"Running PyFrag Driver. Run file: {run_file_path.name}")

    if log_level != logging.DEBUG:
        input_file_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
