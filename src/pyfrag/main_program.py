import argparse
import logging
from pathlib import Path
from typing import Union

from pyfrag import initialize_pyfrag_program
from pyfrag.process_input.process_ams import process_ams_input
from pyfrag.process_input.process_coordfile import extract_molecules_from_coord_file, split_trajectory_into_fragment_molecules
from pyfrag.read_input.parse_input_file import extract_section_blocks_from_file_content
from pyfrag.read_input.read_pyfrag_section import extract_pyfrag_section


def _read_run_file_content(file_path: Union[str, Path]) -> str:
    """Reads the content of a given file. May raises errors when the file does not exist."""
    file = Path(file_path)

    if not file.exists():
        raise FileNotFoundError(f"File {file_path} does not exist. Please make sure the file exists.")

    with open(file, "r") as f:
        return f.read()


def main():
    parser = argparse.ArgumentParser(description="Parse PyFrag input file.")
    parser.add_argument("-f", "--file", help="Path to the run file.", required=True)
    parser.add_argument("-v", "--verbose", help="Set the logger to debug mode.", action="store_true")
    args = parser.parse_args()

    file_path = Path(args.file).resolve()
    log_level = logging.getLevelName("DEBUG") if args.verbose else logging.INFO

    initialize_pyfrag_program(log_level=log_level)
    file_content = _read_run_file_content(file_path)
    file_blocks = extract_section_blocks_from_file_content(file_content)

    if file_blocks.PYFRAG is not None:
        pyfrag_keys = extract_pyfrag_section(file_blocks.PYFRAG)
        print(pyfrag_keys)

    # ==========================================================================
    # ================== Process the data from the input file ==================
    # ==========================================================================

    # ********************* Reading the coordinate file ***********************
    coordinates_data = extract_molecules_from_coord_file(pyfrag_keys["coordfile"])
    trajectories = split_trajectory_into_fragment_molecules(coordinates_data, pyfrag_keys["fragments"])

    # ********************* Converting AMS input to plams settings ***********************
    setting_blocks, calc_type = process_ams_input(file_blocks)


if __name__ == "__main__":
    main()
