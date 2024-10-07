from pathlib import Path

import numpy as np
import pytest
from pyfrag.process_input.process_coordfile import _create_molecule_from_ams_rkf, _create_molecule_from_xyz, _find_files_with_wildcard_option
from scm.plams import Molecule


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


def test_get_molecule_from_pes_rkf(pes_rkf_file_path):
    pes_mols = _create_molecule_from_ams_rkf(pes_rkf_file_path)
    assert len(pes_mols) == 15
    assert isinstance(pes_mols[0], Molecule)


def test_get_molecule_from_xyz(xyz_file_path):
    mols = _create_molecule_from_xyz(xyz_file_path)
    assert len(mols) == 33
    assert isinstance(mols[0], Molecule)


def test_get_molecule_from_amv(amv_file_path):
    mols = _create_molecule_from_xyz(amv_file_path)
    assert len(mols) == 33
    assert isinstance(mols[0], Molecule)


def test_get_molecule_from_irc_rkf(irc_rkf_file_path):
    irc_mols = _create_molecule_from_ams_rkf(irc_rkf_file_path)
    assert len(irc_mols) == 33
    assert isinstance(irc_mols[0], Molecule)


def test_same_molecule_from_irc_xyz_amv(irc_rkf_file_path, xyz_file_path, amv_file_path):
    irc_mols = _create_molecule_from_ams_rkf(irc_rkf_file_path)
    xyz_mols = _create_molecule_from_xyz(xyz_file_path)
    amv_mols = _create_molecule_from_xyz(amv_file_path)

    assert len(irc_mols) == len(xyz_mols) == len(amv_mols), "Lengths are not equal"

    assert np.allclose(irc_mols[0].as_array(), xyz_mols[0].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"
    assert np.allclose(irc_mols[0].as_array(), amv_mols[0].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"
    assert np.allclose(xyz_mols[0].as_array(), amv_mols[0].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"

    assert np.allclose(irc_mols[31].as_array(), xyz_mols[31].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"
    assert np.allclose(irc_mols[31].as_array(), amv_mols[31].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"
    assert np.allclose(xyz_mols[31].as_array(), amv_mols[31].as_array(), rtol=1e-5, atol=1e-8), "Arrays are not close enough"


def test_find_files_with_wildcard_option(mock_filesystem):
    # Test with wildcard
    coord_files = [mock_filesystem / "*.xyz"]
    result = _find_files_with_wildcard_option(coord_files)
    expected = sorted([mock_filesystem / "file1.xyz", mock_filesystem / "file2.xyz"])
    assert sorted(result) == expected

    # Test without wildcard
    coord_files = [mock_filesystem / "file3.amv"]
    result = _find_files_with_wildcard_option(coord_files)
    expected = [mock_filesystem / "file3.amv"]
    assert result == expected

    # Test with wildcard in subdirectory
    coord_files = [mock_filesystem / "subdir" / "*.xyz"]
    result = _find_files_with_wildcard_option(coord_files)
    expected = [mock_filesystem / "subdir" / "file4.xyz"]
    assert result == expected

    # Test with no matching files
    coord_files = [mock_filesystem / "*.nonexistent"]
    result = _find_files_with_wildcard_option(coord_files)
    expected = []
    assert result == expected

    # Test with multiple patterns
    coord_files = [mock_filesystem / "*.xyz", mock_filesystem / "*.amv"]
    result = _find_files_with_wildcard_option(coord_files)
    expected = sorted([mock_filesystem / "file1.xyz", mock_filesystem / "file2.xyz", mock_filesystem / "file3.amv"])
    assert sorted(result) == expected
