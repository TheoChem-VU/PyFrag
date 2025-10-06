import pytest

import pyfrag.executables.adf.input as input_reader
from pyfrag.executables.adf.errors import PyFragSectionInputError


def test_line_with_commas_to_spaces():
    """Note: commas are not supported in the input and are hence replaced with spaces. All the other functions below assume that the input is space-separated and hence are only tested for that."""
    pyfrag_section = """
vdd 0,1,2,3,4
"""
    expected_output = {"vdd": [[0, 1, 2, 3, 4]]}

    assert input_reader.extract_pyfrag_section(pyfrag_section) == expected_output


def test_check_commented_line_valid():
    line1 = "bondlength 1 2 1.5  # With a comment"
    line2 = "bondlength 1 2 1.5  ! With another comment"
    line3 = "bondlength 1 2 1.5  ; With yet another comment!"
    input_key = "bondlength"
    limits = (3, 4)
    expected_output = ["bondlength", "1", "2", "1.5"]
    assert input_reader._check_line_length(line1, input_key, limits) == expected_output
    assert input_reader._check_line_length(line2, input_key, limits) == expected_output
    assert input_reader._check_line_length(line3, input_key, limits) == expected_output


def test_check_line_with_spaces_valid():
    line1 = "    bondlength 1 2 1.5      "
    line2 = "bondlength      1 2 1.5   "
    line3 = "bondlength 1 2        1.5              "
    input_key = "bondlength"
    limits = (3, 4)
    expected_output = ["bondlength", "1", "2", "1.5"]
    assert input_reader._check_line_length(line1, input_key, limits) == expected_output
    assert input_reader._check_line_length(line2, input_key, limits) == expected_output
    assert input_reader._check_line_length(line3, input_key, limits) == expected_output


def test_check_only_comment_line_valid():
    line1 = "# This whole line is a comment "
    line2 = "! This whole line is a comment "
    line3 = "; This whole line is a comment "
    input_key = "bondlength"
    limits = (0, 0)
    expected_output = []
    assert input_reader._check_line_length(line1, input_key, limits) == expected_output
    assert input_reader._check_line_length(line2, input_key, limits) == expected_output
    assert input_reader._check_line_length(line3, input_key, limits) == expected_output


def test_check_line_length_valid():
    line = "bondlength 1 2 1.5"
    input_key = "bondlength"
    limits = (3, 4)
    expected_output = ["bondlength", "1", "2", "1.5"]
    assert input_reader._check_line_length(line, input_key, limits) == expected_output


def test_check_line_length_invalid():
    line = "bondlength 1 2 3 4"
    input_key = "bondlength"
    limits = (3, 4)
    with pytest.raises(PyFragSectionInputError, match="bondlength is not valid. Length of the bondlength is not correct. Make sure to specify the correct format"):
        input_reader._check_line_length(line, input_key, limits)


def test_read_bondlength_line_no_bondlength():
    line = "bondlength 1 2"
    expected_output = (1, 2, 0.0)
    assert input_reader._read_bondlength_line(line) == expected_output


def test_read_bondlength_line_with_bondlength():
    line = "bondlength 1 2 1.5"
    expected_output = (1, 2, 1.5)
    assert input_reader._read_bondlength_line(line) == expected_output


def test_read_bondangle_line_no_angle():
    line = "angle 1 2 3"
    expected_output = (1, 2, 3, 0.0)
    assert input_reader._read_bondangle_line(line) == expected_output


def test_read_bondangle_line_with_angle():
    line = "angle 1 2 3 120.0"
    expected_output = (1, 2, 3, 120.0)
    assert input_reader._read_bondangle_line(line) == expected_output


def test_read_dihedral_angle_no_angle():
    line = "dihedral 1 2 3 4"
    expected_output = (1, 2, 3, 4, 0.0)
    assert input_reader._read_dihedral_angle(line) == expected_output


def test_read_dihedral_angle_with_angle():
    line = "dihedral 1 2 3 4 45.0"
    expected_output = (1, 2, 3, 4, 45.0)
    assert input_reader._read_dihedral_angle(line) == expected_output


def test_read_overlap_line_no_irrep():
    line = "overlap frag1 HOMO frag2 LUMO"
    expected_output = ("frag1", "HOMO", "frag2", "LUMO")
    assert input_reader._read_overlap_line(line) == expected_output


def test_read_overlap_line_with_irrep():
    line = "overlap S frag1 5 AA frag2 4"
    expected_output = ("S", "frag1", "5", "AA", "frag2", "4")
    assert input_reader._read_overlap_line(line) == expected_output


def test_read_population_line_two_fragments_two_MOs():
    line = "population frag1 HOMO"
    expected_output = ("frag1", "HOMO")
    assert input_reader._read_population_line(line) == expected_output


def test_read_population_line_two_fragments_one_MO():
    line = "population frag2 HOMO-1"
    expected_output = ("frag2", "HOMO-1")
    assert input_reader._read_population_line(line) == expected_output


def test_read_population_line_with_irrep():
    line = "population AA frag2 5"
    expected_output = ("AA", "frag2", "5")
    assert input_reader._read_population_line(line) == expected_output


def test_read_orbitalenergy_line_two_fragments_two_MOs():
    line = "orbitalenergy frag1 HOMO"
    expected_output = ("frag1", "HOMO")
    assert input_reader._read_orbitalenergy_line(line) == expected_output


def test_read_orbitalenergy_line_two_fragments_one_MO():
    line = "orbitalenergy frag1 HOMO-2"
    expected_output = ("frag1", "HOMO-2")
    assert input_reader._read_orbitalenergy_line(line) == expected_output


def test_read_orbitalenergy_line_with_irrep():
    line = "ORBITALENERGY AA frag2 5"
    expected_output = ("AA", "frag2", "5")
    assert input_reader._read_orbitalenergy_line(line) == expected_output


def test_read_vdd_line():
    line = "vdd 1 2 3 4 #A comment"
    expected_output = [1, 2, 3, 4]
    assert input_reader._read_vdd_line(line) == expected_output


def test_read_irrep_line():
    line = "  irrepOI        AA #Hello"
    expected_output = ["AA"]
    assert input_reader._read_irrep_line(line) == expected_output


def test_read_strain_line_valid():
    line1 = "strain 0.5"
    line2 = "strain -0.5"
    expected_output1 = 0.5
    expected_output2 = -0.5
    assert input_reader._read_strain_line(line1) == expected_output1
    assert input_reader._read_strain_line(line2) == expected_output2


def test_read_strain_line_invalid():
    line = "strain not_a_number"
    with pytest.raises(PyFragSectionInputError) as excinfo:
        input_reader._read_strain_line(line)
    assert "strain is not valid. Make sure to specify the strain value correctly with a float number" in str(excinfo.value)


def test_read_fragment_indices_line_valid():
    line1 = "fragment 1 2 3 4"
    line2 = "fragment 3-8"
    line3 = "fragment 1 2 3 4 5-8"
    expected_output1 = [1, 2, 3, 4]
    expected_output2 = [3, 4, 5, 6, 7, 8]
    expected_output3 = [1, 2, 3, 4, 5, 6, 7, 8]
    assert input_reader._read_fragment_indices_line(line1) == expected_output1
    assert input_reader._read_fragment_indices_line(line2) == expected_output2
    assert input_reader._read_fragment_indices_line(line3) == expected_output3


def test_read_fragment_indices_line_invalid():
    line1 = "fragment 0"
    line2 = "fragment"
    with pytest.raises(PyFragSectionInputError, match="fragment is not valid. The atom indices should start with 1."):
        input_reader._read_fragment_indices_line(line1)
    with pytest.raises(PyFragSectionInputError, match="Length of the fragment is not correct. Make sure to specify the correct format"):
        input_reader._read_fragment_indices_line(line2)
