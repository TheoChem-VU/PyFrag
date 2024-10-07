import logging
from pathlib import Path
from typing import Callable, Dict, List, Set, Union

import numpy as np
from numpy.typing import NDArray
from scm.plams import Atom, KFHistory, KFReader, Molecule, PeriodicTable, Units

import pyfrag
from pyfrag.errors import PyFragCoordFileError

logger = logging.getLogger(name="Coordfile Processor")
# =============================================================================
# AMV / XYZ file reading  =====================================================
# =============================================================================


def _find_files_with_wildcard_option(coord_files: List[Path]) -> List[Path]:
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
# AMV / XYZ file reading  =====================================================
# =============================================================================


def _create_molecule_from_xyz(xyz_file: Path) -> List[Molecule]:
    """
    Read a molecule or multiple molecules from a .xyz file.
    """
    mols = []

    current_molecule = Molecule()

    lines = xyz_file.read_text().splitlines()

    for line in lines:
        line = line.strip()
        if not line:  # Blank line indicates separation between molecules
            if current_molecule.atoms:  # Check if the molecule has any atoms
                mols.append(current_molecule)
                current_molecule = Molecule()  # Start a new molecule
        else:
            # Assuming the line format is "Symbol X Y Z"
            parts = line.split()

            if parts[0] not in PeriodicTable.symtonum:  # Check if the first part is a valid element symbol
                continue

            if len(parts) == 4:  # Check if the line has 4 parts
                symbol, x, y, z = parts
                atom = Atom(symbol=symbol, coords=(float(x), float(y), float(z)))
                current_molecule.add_atom(atom)

    # After reading the file, check if the last molecule has atoms and hasn't been added yet
    if current_molecule.atoms:
        mols.append(current_molecule)

    return mols


def _create_molecule_from_amv(amv_file: Path) -> List[Molecule]:
    """
    Read a molecule or multiple molecules from an .amv file. Structure is very similar to .xyz file, but with an additional header line for each molecule.
    """
    mols = _create_molecule_from_xyz(amv_file)
    return mols


# =============================================================================
# RKF file reading PESScans and IRCs ==========================================
# =============================================================================


def _create_molecule_from_coords_and_symbols(atom_symbols: List[str], atom_coords: Union[List[List[float]], NDArray]) -> Molecule:
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

    for i, coord in enumerate(history.iter("Coords"), start=1):  # type: ignore  # Results does not have proper typing
        if i not in converged_indices:
            continue

        coord = np.array(coord, dtype=float).reshape(-1, 3) * Units.conversion_ratio("bohr", "angstrom")

        mols.append(_create_molecule_from_coords_and_symbols(atom_symbols, coord))

    return mols


def _extract_molecules_from_irc_path(history: KFHistory, atom_symbols: List[str]) -> List[Molecule]:
    """Extract the molecules from an IRC calculation"""
    mols_forward: List[Molecule] = []
    mols_backward: List[Molecule] = []

    for coord, direction, is_converged in zip(history.iter("Coords"), history.iter("IRCDirection"), history.iter("Converged")):  # type: ignore  # Results does not have proper typing
        if not is_converged:
            continue

        coord = np.array(coord, dtype=float).reshape(-1, 3) * Units.conversion_ratio("bohr", "angstrom")
        mol = _create_molecule_from_coords_and_symbols(atom_symbols, coord)
        if direction == 1:  # Forward
            mols_forward.append(mol)
        elif direction == 2:  # Backward
            mols_backward.append(mol)
        else:
            raise ValueError(f"Unknown IRC direction: {direction}")
    return mols_backward[::-1] + mols_forward


def _create_molecule_from_ams_rkf(rkf_file: Path) -> List[Molecule]:
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


extension_func_mapping: Dict[str, Callable[[Path], List[Molecule]]] = {
    "xyz": _create_molecule_from_xyz,
    "amv": _create_molecule_from_amv,
    "rkf": _create_molecule_from_ams_rkf,
}


def extract_molecules_from_coord_file(coord_files: List[Path]) -> List[Molecule]:
    """
    Extract molecules from one or multiple coordinate files.
    Supported file types are .xyz, .amv, .rkf
    """

    coord_files = _find_files_with_wildcard_option(coord_files)

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
# Fragmentation of molecules  =================================================
# =============================================================================


def split_trajectory_into_fragment_molecules(mols: List[Molecule], frag_indices: List[List[int]]) -> List[List[Molecule]]:
    """
    Given a list of molecules, split each molecule into fragments (based on the entries = number of fragments in the nested frag_indices list).

    Returns a list with length n_frags + 1: [[complex_trajectory], [frag1_trajectory], [frag2_trajectory], ...]
    where the first entry is the trajectory of the complex molecule and the following entries are the trajectories of the fragments
    """
    all_specified_indices = _collect_specified_indices(frag_indices=frag_indices)
    frag_indices = _handle_special_indices_case(n_atoms=len(mols[0]), frag_indices=frag_indices, all_specified_indices=all_specified_indices)

    logger.debug(f"Fragment indices: {frag_indices}")

    trajectories = _create_trajectory_lists(mols=mols, frag_indices=frag_indices)
    return trajectories


def _collect_specified_indices(frag_indices: List[List[int]]) -> Set[int]:
    """
    Collect all explicitly specified indices from the frag_indices list, excluding -1.
    """
    all_specified_indices = set()
    for indices in frag_indices:
        if indices != [-1]:
            all_specified_indices.update(indices)
    return all_specified_indices


def _handle_special_indices_case(n_atoms: int, frag_indices: List[List[int]], all_specified_indices: Set[int]) -> List[List[int]]:
    """
    Handle the special case where frag_indices contains [-1]. This indicates that all indices that have not been specified should be used.
    An example is [[-1], [1, 2, 3]] -> [[4, 5, 6], [1, 2, 3]]
    """
    for i, indices in enumerate(frag_indices):
        if indices == [-1]:
            all_indices = set(range(1, n_atoms + 1))
            remaining_indices = all_indices - all_specified_indices
            frag_indices[i] = list(remaining_indices)
    return frag_indices


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


def main():
    coord_file_dir = Path(pyfrag.__file__).parent.parent.parent.resolve() / "tests" / "fixtures" / "coord_files"
    fragment_indices = [[-1], [3, 5, 6]]

    files = [
        ("xyz", coord_file_dir / "ams.xyz"),
        ("amv", coord_file_dir / "ams.amv"),
        ("rkf", coord_file_dir / "pes.ams.rkf"),
        ("rkf", coord_file_dir / "irc.ams.rkf"),
    ]

    for pair in files:
        mol = Molecule()
        ext, file = pair
        mol = extension_func_mapping[ext](file)
        res = split_trajectory_into_fragment_molecules(mol, fragment_indices)

        print(res)


if __name__ == "__main__":
    main()
