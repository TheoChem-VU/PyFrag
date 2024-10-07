import argparse
import logging
from pathlib import Path
from typing import Union

from pyfrag import initialize_pyfrag_program
from pyfrag.jobs.basejob import RestrictedPyFragJob
from pyfrag.jobs.fragment import create_fragments_and_complex
from pyfrag.process_input.process_coordfile import extract_molecules_from_coord_file, split_trajectory_into_fragment_molecules
from pyfrag.read_input.parse_input_file import extract_section_blocks_from_file_content
from pyfrag.read_input.read_pyfrag_section import extract_pyfrag_section


def read_input_file_content(file_path: Union[str, Path]) -> str:
    """Reads the content of a given file. May raises errors when the file does not exist."""
    file = Path(file_path)

    if not file.exists():
        raise FileNotFoundError(f"File {file_path} does not exist. Please make sure the file exists.")

    with open(file, "r") as f:
        lines = f.read()
    return lines


class TestSystems:
    symmetry_test: Path = Path(__file__).parent.parent.parent / "example" / "restricted" / "developer_test_symmetry" / "symmetry_job.in"
    restricted_test: Path = Path(__file__).parent.parent.parent / "example" / "restricted" / "regular_job" / "new_ams_job.in"
    irreps_test: Path = Path(__file__).parent.parent.parent / "example" / "restricted" / "developer_test_irreps" / "irrep_job.in"
    charged_test: Path = Path(__file__).parent.parent.parent / "example" / "restricted" / "developer_test_charged" / "PES1Dist57Angle80.in"
    unrestricted_test: Path = Path(__file__).parent.parent.parent / "example" / "open-shell" / "new_ams_methane.in"


def main():
    # parser = argparse.ArgumentParser(description="Parse PyFrag input file.")
    # parser.add_argument("-f", "--file", help="Path to the run file.", required=True)
    # parser.add_argument("-v", "--verbose", help="Enable the debug option.", required=False, action="store_true")
    # args = parser.parse_args()

    # args = argparse.Namespace(file=str(Path(__file__).parent.parent.parent / "example" / "restricted" / "developer_test_symmetry" / "symmetry_job.in"), verbose=True)
    args = argparse.Namespace(file=str(TestSystems.symmetry_test), verbose=False)

    file_path = Path(args.file).resolve()
    log_level = logging.DEBUG if args.verbose else logging.INFO

    initialize_pyfrag_program(log_level=log_level)

    # ==========================================================================
    # ================== Reading the input file content ========================
    # ==========================================================================

    file_content = read_input_file_content(file_path)
    file_blocks = extract_section_blocks_from_file_content(file_content)
    file_blocks.validate()  # Checks if the input file has the necessary sections such as PyFrag, ADF, and AMS

    pyfrag_keys = extract_pyfrag_section(file_blocks.PYFRAG, file_path.parent)

    # ==========================================================================
    # ================== Create the Fragment and Complex =======================
    # ==========================================================================

    coordinates_data = extract_molecules_from_coord_file(pyfrag_keys.coordfile)
    mol_trajectories = split_trajectory_into_fragment_molecules(coordinates_data, pyfrag_keys.fragment_indices)

    complex, fragments = create_fragments_and_complex(file_blocks=file_blocks, mol_trajectories=mol_trajectories, complex_name="complex", frag_names=["frag1", "frag2"], pyfrag_keys=pyfrag_keys)

    # ==========================================================================
    # ================== Create and Run the PyFragJob ==========================
    # ==========================================================================

    exit()
    with RestrictedPyFragJob(job_dir=file_path.parent, complex=complex, fragments=fragments) as job:
        job.run()

    # ==========================================================================
    # ================== Analyse the results ===================================
    # ==========================================================================


if __name__ == "__main__":
    main()
