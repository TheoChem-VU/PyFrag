from pathlib import Path
from typing import Callable, Dict, List, Union

import numpy as np
from numpy.typing import NDArray
from scm.plams import Atom, KFHistory, KFReader, Molecule, PeriodicTable, Units

import pyfrag

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


def main():
    coord_file_dir = Path(pyfrag.__file__).parent.parent.parent.resolve() / "tests" / "fixtures" / "coord_files"

    extensions = ["xyz", "rkf", "amv"]
    files = {ext: coord_file_dir / f"ams.{ext}" for ext in extensions}

    for ext, file in files.items():
        mol = Molecule()
        mol = extension_func_mapping[ext](file)
        print(ext, len(mol), mol[31])


if __name__ == "__main__":
    main()
