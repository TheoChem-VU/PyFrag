from pathlib import Path

import numpy as np
import pytest
from pyfrag.process_input.read_coordfile import _create_molecule_from_ams_rkf, _create_molecule_from_xyz
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
