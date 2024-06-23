from pathlib import Path
from typing import Callable, Dict

from scm.plams import Molecule

import pyfrag


def _create_molecule_from_xyz(mol: Molecule, xyz_file: Path):
    """
    Read a molecule or multiple molecules from a .xyz file.
    """
    mol.read(xyz_file, "xyz")  # type: ignore  # plams does not have static typing
    return mol


def _create_molecule_from_ams_rkf(mol: Molecule, rkf_file: Path):
    """
    Read a molecule or multiple molecules from an AMS .rkf file.
    """
    mol = mol.read(rkf_file, "rkf")  # type: ignore  # plams does not have static typing
    return mol


def _create_molecule_from_amv(mol: Molecule, amv_file: Path):
    """
    Read a molecule or multiple molecules from an .amv file
    """
    mol = mol.read(amv_file, "xyz")  # type: ignore  # plams does not have static typing
    return mol


extension_func_mapping: Dict[str, Callable[[Molecule, Path], Molecule]] = {
    "xyz": _create_molecule_from_xyz,
    "amv": _create_molecule_from_amv,
    "rkf": _create_molecule_from_ams_rkf,
}


def main():
    coord_file_dir = Path(pyfrag.__file__).parent.parent.parent.resolve() / "tests" / "fixtures" / "coord_files"
    print(coord_file_dir)

    extensions = ["xyz", "rkf", "amv"]
    files = {ext: coord_file_dir / f"ams.{ext}" for ext in extensions}

    for ext, file in files.items():
        mol = Molecule()
        mol = extension_func_mapping[ext](mol, file)
        print(ext, mol)


if __name__ == "__main__":
    main()
