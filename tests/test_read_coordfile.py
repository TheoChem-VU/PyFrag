from pathlib import Path

import pytest
from pyfrag.process_input.read_coordfile import _create_molecule_from_ams_rkf, _create_molecule_from_amv, _create_molecule_from_xyz
from scm.plams import Molecule


@pytest.fixture(params=["xyz", "amv", "rkf"])
def coord_file_path(request):
    current_file_path = Path(__file__).parent.resolve()
    return current_file_path / "fixtures" / "coord_files" / f"ams.{request.param}"


def test_get_molecule_from_xyz(coord_file_path):
    if coord_file_path.suffix == ".xyz":
        mol = Molecule()
        result = _create_molecule_from_xyz(mol, coord_file_path)
        # Replace the following line with actual validation for your Molecule object
        assert result is not None, "Failed to read .xyz file"


def test_create_molecule_from_ams_rkf(coord_file_path):
    if coord_file_path.suffix == ".rkf":
        mol = Molecule()
        result = _create_molecule_from_ams_rkf(mol, coord_file_path)
        # Replace the following line with actual validation for your Molecule object
        assert result is not None, "Failed to read .rkf file"


def test_create_molecule_from_amv(coord_file_path):
    if coord_file_path.suffix == ".amv":
        mol = Molecule()
        result = _create_molecule_from_amv(mol, coord_file_path)
        # Replace the following line with actual validation for your Molecule object
        assert result is not None, "Failed to read .rkf file"
