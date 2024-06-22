import argparse
from pathlib import Path
from typing import Union

from pyfrag import initialize_pyfrag_program
from pyfrag.input.parse_input_file import extract_section_blocks_from_file_content
from pyfrag.input.read_pyfrag_section import extract_pyfrag_section


def _read_run_file_content(file_path: Union[str, Path]) -> str:
    """Reads the content of a given file. May raises errors when the file does not exist."""
    file = Path(file_path)

    if not file.exists():
        raise FileNotFoundError(f"File {file_path} does not exist. Please make sure the file exists.")

    with open(file, "r") as f:
        return f.read()


def main():
    initialize_pyfrag_program()
    parser = argparse.ArgumentParser(description="Parse PyFrag input file.")
    parser.add_argument("-f", "--file", help="Path to the run file.", required=True)
    args = parser.parse_args()

    file_content = _read_run_file_content(args.file)
    file_blocks = extract_section_blocks_from_file_content(file_content)

    if file_blocks.PYFRAG is not None:
        pyfrag_keys = extract_pyfrag_section(file_blocks.PYFRAG)
        print(pyfrag_keys)


if __name__ == "__main__":
    main()
