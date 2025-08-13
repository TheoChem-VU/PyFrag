from pathlib import Path

import numpy as np
import pytest
from scm.plams import Atom, Molecule

from pyfrag.executables.adf import mol_handling


# Fixture to generate a single Molecule object
@pytest.fixture
def molecule_fixture():
    mol = Molecule()
    atoms = ["H", "O", "N", "C", "P", "S", "F"]
    for i, symbol in enumerate(atoms, start=1):
        mol.add_atom(Atom(symbol=symbol, coords=[i * 0.1, i * 0.1, i * 0.1]))
    return mol


# Fixture to generate a list of Molecule objects for a trajectory using the molecule_fixture
@pytest.fixture
def trajectory_fixture(molecule_fixture: Molecule):
    trajectory = [molecule_fixture.copy() for _ in range(5)]
    return trajectory


@pytest.fixture()
def pes_rkf_file_path():
    current_file_path = Path(__file__).parent.resolve()
    return current_file_path / "fixtures" / "coord_files" / "pes.ams.rkf"


@pytest.fixture()
def xyz_file_path():
    current_file_path = Path(__file__).parent.resolve()
    return current_file_path / "fixtures" / "coord_files" / "ams.xyz"


@pytest.fixture()
def amv_file_path():
    current_file_path = Path(__file__).parent.resolve()
    return current_file_path / "fixtures" / "coord_files" / "ams.amv"


@pytest.fixture()
def irc_rkf_file_path():
    current_file_path = Path(__file__).parent.resolve()
    return current_file_path / "fixtures" / "coord_files" / "irc.ams.rkf"


@pytest.fixture()
def mock_filesystem(tmp_path):
    # Create a temporary directory structure with some files
    (tmp_path / "file1.xyz").touch()
    (tmp_path / "file2.xyz").touch()
    (tmp_path / "file3.amv").touch()
    (tmp_path / "subdir").mkdir()
    (tmp_path / "subdir" / "file4.xyz").touch()
    return tmp_path


# -----------------------------------------------------------------------------
# Unit tests regarding loading molecules and defining fragments
# -----------------------------------------------------------------------------


def test_is_header_line():
    assert mol_handling.is_header_line("Geometry 1, Name: O_SC1_E_c1, Energy: -4528.1 Ha")
    assert not mol_handling.is_header_line("C 0.0 0.0 0.0")
    assert mol_handling.is_header_line("")
    assert not mol_handling.is_header_line(["C", "0.0", "0.0", "0.0"])


def test_load_molecule_from_str():
    block = """C 0.0 0.0 0.0\nH 0.0 0.0 1.0\n"""
    mol = mol_handling.load_molecule_from_str(block)
    assert len(mol.atoms) == 2
    assert mol.atoms[0].symbol == "C"
    assert np.allclose(mol.atoms[1].coords, [0.0, 0.0, 1.0])


def test_expand_indices():
    # _expand_indices is private, but we can test via update_fragment_indices
    indices = [[1, "2-4"], [5]]
    expanded = mol_handling.update_fragment_indices(indices)
    assert expanded == [[1, 2, 3, 4], [5]]


def test_update_fragment_indices_uniqueness():
    indices = [[1, 2], [3, 4]]
    expanded = mol_handling.update_fragment_indices(indices)
    assert expanded == [[1, 2], [3, 4]]
    # Overlapping indices should raise
    with pytest.raises(Exception):
        mol_handling.update_fragment_indices([[1, 2], [2, 3]])


def test_validate_unique_indices(molecule_fixture):
    # Should not raise
    mol_handling.validate_unique_indices([[1, 2], [3, 4]], molecule=molecule_fixture)
    # Out of range should raise
    with pytest.raises(Exception):
        mol_handling.validate_unique_indices([[1, 8]], molecule=molecule_fixture)


def test_split_trajectory_into_fragment_molecules(trajectory_fixture):
    # 7 atoms in molecule_fixture, so indices 1,2,3 are valid
    frag_indices = [[1, 2], [3]]
    result = mol_handling.split_trajectory_into_fragment_molecules(trajectory_fixture, frag_indices)
    # result[0] is the complex trajectory, result[1] and result[2] are fragments
    assert len(result) == 3
    assert all(isinstance(m, type(trajectory_fixture[0])) for m in result[0])
    assert all(len(frag.atoms) == 2 for frag in result[1])
    assert all(len(frag.atoms) == 1 for frag in result[2])


# -----------------------------------------------------------------------------
# Unit tests regarding loading molecules from various file types
# -----------------------------------------------------------------------------


def test_get_molecule_from_pes_rkf(pes_rkf_file_path):
    pes_mols = mol_handling.create_molecules_from_rkf_file(pes_rkf_file_path)
    assert len(pes_mols) == 15
    assert isinstance(pes_mols[0], Molecule)


def test_get_molecule_from_xyz(xyz_file_path):
    mols = mol_handling.create_molecules_from_xyz_file(xyz_file_path)
    assert len(mols) == 33
    assert isinstance(mols[0], Molecule)


def test_get_molecule_from_amv(amv_file_path):
    mols = mol_handling.create_molecules_from_amv_file(amv_file_path)
    assert len(mols) == 33
    assert isinstance(mols[0], Molecule)


def test_get_molecule_from_irc_rkf(irc_rkf_file_path):
    irc_mols = mol_handling.create_molecules_from_rkf_file(irc_rkf_file_path)
    assert len(irc_mols) == 33
    assert isinstance(irc_mols[0], Molecule)


def test_same_molecule_from_irc_xyz_amv(irc_rkf_file_path, xyz_file_path, amv_file_path):
    irc_mols = mol_handling.create_molecules_from_rkf_file(irc_rkf_file_path)
    xyz_mols = mol_handling.create_molecules_from_xyz_file(xyz_file_path)
    amv_mols = mol_handling.create_molecules_from_amv_file(amv_file_path)

    assert len(irc_mols) == len(xyz_mols) == len(amv_mols), "Lengths are not equal"

    assert np.allclose(irc_mols[0].as_array(), xyz_mols[0].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"
    assert np.allclose(irc_mols[0].as_array(), amv_mols[0].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"
    assert np.allclose(xyz_mols[0].as_array(), amv_mols[0].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"

    assert np.allclose(irc_mols[31].as_array(), xyz_mols[31].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"
    assert np.allclose(irc_mols[31].as_array(), amv_mols[31].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"
    assert np.allclose(xyz_mols[31].as_array(), amv_mols[31].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"


def testfind_files_with_wildcard_option(mock_filesystem):
    # Test with wildcard
    coord_files = [mock_filesystem / "*.xyz"]
    result = mol_handling.find_files_with_wildcard_option(coord_files)
    expected = sorted([mock_filesystem / "file1.xyz", mock_filesystem / "file2.xyz"])
    assert sorted(result) == expected

    # Test without wildcard
    coord_files = [mock_filesystem / "file3.amv"]
    result = mol_handling.find_files_with_wildcard_option(coord_files)
    expected = [mock_filesystem / "file3.amv"]
    assert result == expected

    # Test with wildcard in subdirectory
    coord_files = [mock_filesystem / "subdir" / "*.xyz"]
    result = mol_handling.find_files_with_wildcard_option(coord_files)
    expected = [mock_filesystem / "subdir" / "file4.xyz"]
    assert result == expected

    # Test with no matching files
    coord_files = [mock_filesystem / "*.nonexistent"]
    result = result = mol_handling.find_files_with_wildcard_option(coord_files)
    expected = []
    assert result == expected

    # Test with multiple patterns
    coord_files = [mock_filesystem / "*.xyz", mock_filesystem / "*.amv"]
    result = mol_handling.find_files_with_wildcard_option(coord_files)
    expected = sorted([mock_filesystem / "file1.xyz", mock_filesystem / "file2.xyz", mock_filesystem / "file3.amv"])
    assert sorted(result) == expected
