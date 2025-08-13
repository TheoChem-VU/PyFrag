import logging
import pathlib as pl
import re
from typing import TYPE_CHECKING, Callable, Dict, List, Sequence, Union

import numpy as np
from scm.plams import AMSJob, Atom, KFHistory, KFReader, Molecule

# Since we use amspython (from the AMS pacakge) as python environment, we cannot install the pyfrag package.
# However, when performing tests with pytest, we install the pyfrag package and then the "." needs to be added.
# Therefore, we include this try ... except block.
try:
    from .constants import BOHR_TO_ANGSTROM
    from .errors import FragmentIndicesError, PyFragCoordFileError
except ImportError:
    from constants import BOHR_TO_ANGSTROM
    from errors import FragmentIndicesError, PyFragCoordFileError


if TYPE_CHECKING:
    from numpy.typing import NDArray


logger = logging.getLogger(name="Coordfile Processing-ADF")
logger.setLevel(logging.DEBUG)


# =============================================================================
# File finding with wildcard option (e.g. *.xyz)
# =============================================================================


def find_files_with_wildcard_option(coord_files: Sequence[pl.Path]) -> List[pl.Path]:
    """
    Find files with the * option ("wildcard") in the path.
    """
    files = []
    for coord_file in coord_files:
        if "*" in coord_file.name:
            files.extend(coord_file.parent.glob(coord_file.name))
        else:
            files.append(coord_file)
    return files


# =============================================================================
# AMV / XYZ file reading
# =============================================================================


def is_header_line(line: Union[Sequence[str], str]) -> bool:
    """Check if the line of the block is a header line such as in .amv files:

    Geometry 1, Name: O_SC1_E_c1, Energy: -4528.147256459182 Ha

    which is a header line for the first molecule in the block.

    Args:
        str_block (list[str]): A list of strings representing the lines in the block.
    Returns:
        bool: True if the first line is a header line, False otherwise.
    """

    if not line:
        return True

    parts = line.split() if isinstance(line, str) else line
    # Check if the first line is in the [str, float, float, float, [*other]] format
    try:
        # Check if the first part is a string (name) and the rest are floats (coordinates)
        (_, _, _, _) = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
        return False
    except (ValueError, IndexError):
        # If conversion to float fails, there is a header line
        return True


def load_molecule_from_str(molecule_block: Union[str, Sequence[str]]) -> Molecule:
    """Load a molecule from a string by adding atoms and coordinates explicitly.

    Args:
        molecule_block (str | list[str]): The molecule block as a string or list of strings.

    Returns:
        Molecule: A PLAMS Molecule object.
    """
    if isinstance(molecule_block, str):
        molecule_block = molecule_block.split("\n")

    # Remove empty lines and strip whitespace
    molecule_block = [line.strip() for line in molecule_block if line.strip()]

    # Create a new Molecule object
    mol = Molecule()

    # Iterate over the lines in the molecule block
    for line in molecule_block:
        # Example line: "O 0.000000 0.000000 0.000000 region=O.adf"
        parts = line.split()
        if len(parts) >= 4:
            atom_symbol: str = parts[0]
            coordinates: list[float] = [float(coord) for coord in parts[1:4]]
            atom = Atom(symbol=atom_symbol, coords=coordinates)
            mol.add_atom(atom)

            # Check if there are additional properties (e.g., region details) which need to be added to the suffix key in the atom.properties Settings
            # Example: "region=O.adf"
            if len(parts) > 4:
                atom.properties.suffix = " ".join(parts[4:])
    return mol


def create_molecules_from_xyz_file(xyz_file: pl.Path) -> List[Molecule]:
    """
    Read a molecule or multiple molecules from a .xyz file.
    """
    file_content = xyz_file.read_text()

    # Extract the blocks which are seperated by empty lines
    blocks = re.split(r"\n\s*\n", file_content.strip())

    molecule_strings: list[list[str]] = []
    for mol_block in blocks:
        # Check if the line is a header line
        lines = mol_block.splitlines()

        # Skip empty lines
        if not lines:
            continue

        if is_header_line(lines[0]):
            # If the first line is a header line, we need to skip it and add the rest of the lines
            lines = lines[1:]

        molecule_strings.append(lines)

    molecules = [load_molecule_from_str(mol_str) for mol_str in molecule_strings]

    return molecules


def create_molecules_from_amv_file(amv_file: pl.Path) -> List[Molecule]:
    """
    Read a molecule or multiple molecules from an .amv file. Structure is very similar to .xyz file, but with an additional header line for each molecule.
    """
    mols = create_molecules_from_xyz_file(amv_file)
    return mols


# =============================================================================
# RKF file reading PESScans and IRCs
# =============================================================================


def _create_molecule_from_coords_and_symbols(atom_symbols: Sequence[str], atom_coords: Union[Sequence[Sequence[float]], "NDArray"]) -> Molecule:
    """
    Create a PLAMS molecule from atom symbols and coordinates.
    """
    mol = Molecule()
    for atom_symbol, atom_coord in zip(atom_symbols, atom_coords):
        mol.add_atom(Atom(symbol=atom_symbol, coords=atom_coord))
    return mol


def _extract_molecules_from_pes_scan(history: KFHistory, atom_symbols: List[str], converged_indices: List[float]) -> List[Molecule]:
    """Extract the molecules from a PES scan calculation"""
    mols: List[Molecule] = []

    for i, coordinates in enumerate(history.iter("Coords"), start=1):  # type: ignore  # Results does not have proper typing
        if i not in converged_indices:
            continue

        coordinates = np.array(coordinates, dtype=float).reshape(-1, 3) * BOHR_TO_ANGSTROM

        mols.append(_create_molecule_from_coords_and_symbols(atom_symbols, coordinates))

    return mols


def _extract_molecules_from_irc_path(history: KFHistory, atom_symbols: List[str]) -> List[Molecule]:
    """Extract the molecules from an IRC calculation"""
    mols_forward: List[Molecule] = []
    mols_backward: List[Molecule] = []

    for coord, direction, is_converged in zip(history.iter("Coords"), history.iter("IRCDirection"), history.iter("Converged")):  # type: ignore  # Results does not have proper typing
        if not is_converged:
            continue

        coord = np.array(coord, dtype=float).reshape(-1, 3) * BOHR_TO_ANGSTROM
        mol = _create_molecule_from_coords_and_symbols(atom_symbols, coord)
        if direction == 1:  # Forward
            mols_forward.append(mol)
        elif direction == 2:  # Backward
            mols_backward.append(mol)
        else:
            raise ValueError(f"Unknown IRC direction: {direction}")
    return mols_backward[::-1] + mols_forward


def create_molecules_from_rkf_file(rkf_file: pl.Path) -> List[Molecule]:
    """
    Read a molecule or multiple molecules from an AMS .rkf file.
    """
    kf_reader = KFReader(rkf_file)
    history_reader = KFHistory(kf_reader, "History")
    atom_symbols: List[str] = kf_reader.read("Molecule", "AtomSymbols").split()  # type: ignore  # Results does not have proper typing

    if ("IRC", "direction") in kf_reader:
        mols = _extract_molecules_from_irc_path(history_reader, atom_symbols)
    elif ("PESScan", "nPoints") in kf_reader:
        converged_points = [index for index in kf_reader.read("PESScan", "HistoryIndices")]  # type: ignore  # Results does not have proper typing
        mols = _extract_molecules_from_pes_scan(history_reader, atom_symbols, converged_points)
    else:
        raise ValueError(f"Unsupported .rkf file: {rkf_file}")

    return mols


extension_func_mapping: Dict[str, Callable[[pl.Path], List[Molecule]]] = {
    "xyz": create_molecules_from_xyz_file,
    "amv": create_molecules_from_amv_file,
    "rkf": create_molecules_from_rkf_file,
}


def extract_molecules_from_coord_file(coord_files: Sequence[pl.Path]) -> List[Molecule]:
    """
    Extract molecules from one or multiple coordinate files.
    Supported file types are .xyz, .amv, .rkf
    """

    coord_files = find_files_with_wildcard_option(coord_files)

    if not all(coord_file.exists() for coord_file in coord_files):
        raise PyFragCoordFileError(f"File(s) not found: {coord_files}")

    extensions = {coord_file.suffix[1:] for coord_file in coord_files}

    if not all(extension in extension_func_mapping for extension in extensions):
        raise PyFragCoordFileError(f"Unsupported file extension detected in: {extensions}. Allowed are {set(extension_func_mapping)}.")

    molecules = []
    for extension, coord_file in zip(extensions, coord_files):
        molecules.extend(extension_func_mapping[extension](coord_file))

    return molecules


# =============================================================================
# Updating and validating fragment indices
# =============================================================================


def update_fragment_indices(fragment_indices: Sequence[Sequence[Union[int, str]]]) -> List[List[int]]:
    """
    Given a list of fragment indices which may contain the range option (e.g., "3-7"),
    this function updates the fragment indices to include all atoms that are not specified in the list.

    Example:
        fragment_indices = [[1, 2], ["3-7"]] will be updated to [[1, 2], [3, 4, 5, 6, 7]]

    Args:
        complex_mol (Molecule): The complex molecule.
        fragment_indices (List[List[Union[int, str]]]): A list of lists containing the indices of the atoms to be used for fragmentation.

    Returns:
        List[List[int]]: The updated list of fragment indices.
    """
    specified_indices = set()
    expanded_indices: List[List[int]] = []

    # Expand all fragment indices
    for indices in fragment_indices:
        expanded = _expand_indices(indices)
        expanded_indices.append(expanded)
        specified_indices.update(expanded)

    # Validate uniqueness of indices
    validate_unique_indices(expanded_indices)

    return expanded_indices


def _expand_indices(indices: Sequence[Union[int, str]]) -> List[int]:
    """
    Expand a list of indices, handling single integers and range options (e.g., "3-7").
    """
    expanded = []
    for index in indices:
        index = int(index) if isinstance(index, str) and index.isdigit() else index

        if isinstance(index, int):
            if index > 0:  # Only add positive indices
                expanded.append(index)
        elif isinstance(index, str) and "-" in index:
            start, end = map(int, index.split("-"))
            expanded.extend([i for i in range(start, end + 1) if i > 0])  # Only add positive indices
    return expanded


def validate_unique_indices(expanded_indices: List[List[int]], molecule: Union[Molecule, None] = None) -> None:
    """
    Ensure that all indices across fragments are unique and match the length of the molecule.
    """
    all_indices = {i for indices in expanded_indices for i in indices}
    if len(all_indices) != sum(len(indices) for indices in expanded_indices):
        raise FragmentIndicesError(f"Fragment indices are not unique. Please check your 'fragment' keys in the PyFrag input file: {expanded_indices}")

    # If the molecule is provided, check if the indices are within the range of the molecule
    if molecule is not None:
        for indices in expanded_indices:
            for index in indices:
                if index < 1 or index > len(molecule.atoms):
                    raise FragmentIndicesError(f"Fragment indices exceed the number of atoms in the molecule. Please check your 'fragment' keys in the PyFrag input file: {expanded_indices}")

        # And check if the number of indices is equal to the number of atoms in the molecule
        for indices in expanded_indices:
            if len(indices) != len(set(indices)):
                raise FragmentIndicesError(f"Fragment indices are not unique. Please check your 'fragment' keys in the PyFrag input file: {expanded_indices}")


# =============================================================================
# Fragmentation of the trajectory into complex and fragments
# =============================================================================


def split_trajectory_into_fragment_molecules(mols: List[Molecule], frag_indices: Sequence[Sequence[Union[int, str]]]) -> List[List[Molecule]]:
    """
    Given a list of molecules, split each molecule into fragments (based on the entries = number of fragments in the nested frag_indices list).

    Returns a list with length n_frags + 1: [[complex_trajectory], [frag1_trajectory], [frag2_trajectory], ...]
    where the first entry is the trajectory of the complex molecule and the following entries are the trajectories of the fragments
    """
    updated_frag_indices = update_fragment_indices(frag_indices)

    validate_unique_indices(updated_frag_indices)
    logger.info(f"Updated and validated fragment indices: {updated_frag_indices}")

    trajectories = _create_trajectory_lists(mols=mols, frag_indices=updated_frag_indices)
    return trajectories


def _create_trajectory_lists(mols: List[Molecule], frag_indices: List[List[int]]) -> List[List[Molecule]]:
    """
    Create the trajectory lists for the complex molecule and fragments.
    """
    trajectories = [mols]
    for indices in frag_indices:
        fragments = []
        for mol in mols:
            fragments.append(mol.copy([mol.atoms[i - 1] for i in indices]))
        trajectories.append(fragments)
    return trajectories


# =============================================================================
# Main function
# =============================================================================


def create_pyfrag_trajectory_from_coord_file(coord_file: Union[Sequence[pl.Path], pl.Path], fragment_indices: Sequence[Sequence[Union[int, str]]]) -> List[List[Molecule]]:
    """
    Create a trajectory from a coordinate file and split it into fragments.

    Args:
        coord_file (pl.Path): The path to the coordinate file. This can be a .xyz, .amv, or .rkf file.
        fragment_indices (List[List[int]]): A list of lists containing the indices of the atoms to be used for fragmentation.

    Returns:
        List[List[Molecule]]: A list of lists containing the trajectories of the complex and fragments. It is nested list in th order of [[complex@point1, complex@point2, ...], [frag1@point1, frag1@point2, ...], [frag2@point1, frag2@point2, ...], ...]

    Note: all the coordinates are always reported to Angstrom.
    """
    # Multiple coordinate files can be passed as a list
    coord_file = [coord_file] if isinstance(coord_file, pl.Path) else coord_file
    coord_file = [pl.Path(file) for file in coord_file]

    # Extract molecules from the coordinate file
    mols = extract_molecules_from_coord_file(coord_file)

    # Split the trajectory into fragments
    trajectories = split_trajectory_into_fragment_molecules(mols, fragment_indices)

    return trajectories


# =============================================================================
# Test function
# =============================================================================


def main_test():
    # coord_file_dir = Path(__file__).parent.parent.parent.resolve() / "tests" / "fixtures" / "coord_files"
    coord_file_dir = pl.Path(__file__).parent / "example"
    fragment_indices = [["2"], ["1", "3-6"]]

    files = [
        # ("xyz", coord_file_dir /  "developer_test_with_regions" / "mol_trajectory.xyz"),
        ("xyz", coord_file_dir / "palladium" / "molecule.xyz"),
        # ("amv", coord_file_dir / "SC-optimizations_optimized_molecules.amv"),
        # ("rkf", coord_file_dir / "pes.ams.rkf"),
        # ("rkf", coord_file_dir / "irc.ams.rkf"),
    ]

    for pair in files:
        mol = Molecule()
        ext, file = pair
        mol = extension_func_mapping[ext](file)
        res = split_trajectory_into_fragment_molecules(mol, fragment_indices)  # type: ignore  # Results does not have proper typing
        job_input_script = AMSJob(name="test", molecule=res[2][0]).get_input()
        print(job_input_script)


if __name__ == "__main__":
    main_test()
