import pytest
from pyfrag.process_input.process_coordfile import _collect_specified_indices, _handle_special_indices_case, split_trajectory_into_fragment_molecules
from scm.plams import Atom, Molecule


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


def test_collect_specified_indices():
    frag_indices = [[1, 2], [-1], [3, 4]]
    expected = {1, 2, 3, 4}
    result = _collect_specified_indices(frag_indices)
    assert result == expected, "The collected indices do not match the expected set."


def test_handle_special_case():
    n_atoms = 20
    frag_indices = [[-1], [1, 2, 3], [16, 17]]
    all_specified_indices = {1, 2, 3, 16, 17}
    expected_output = [[4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20], [1, 2, 3], [16, 17]]
    result = _handle_special_indices_case(n_atoms, frag_indices, all_specified_indices)
    assert result == expected_output, "The function did not return the expected output."


def test_handle_special_case_sorted_order():
    n_atoms = 6
    frag_indices = [[-1], [1, 2, 3]]
    all_specified_indices = {1, 2, 3}
    expected_output = [[4, 5, 6], [1, 2, 3]]

    # Call the function with the test case
    result = _handle_special_indices_case(n_atoms, frag_indices, all_specified_indices)

    # Sort the sublists and the outer list to ensure the order doesn't affect the test
    result_sorted = [sorted(sublist) for sublist in result]
    result_sorted.sort(key=lambda x: x[0])
    expected_output_sorted = [sorted(sublist) for sublist in expected_output]
    expected_output_sorted.sort(key=lambda x: x[0])

    assert result_sorted == expected_output_sorted, "The function did not return the expected output."


def test_split_trajectory_into_fragment_molecules(trajectory_fixture):
    frag_indices = [[1, 2], [3, 4, 5]]
    result = split_trajectory_into_fragment_molecules(trajectory_fixture, frag_indices)
    assert len(result) == 3, "The number of trajectories returned does not match the number of fragments + 1 (complex)."
    assert len(result[0]) == len(trajectory_fixture), "The complex trajectory does not have the same length as the input trajectory."
    assert len(result[1]) == len(trajectory_fixture), "The first fragment trajectory does not have the same length as the input trajectory."
    assert len(result[2]) == len(trajectory_fixture), "The second fragment trajectory does not have the same length as the input trajectory."
