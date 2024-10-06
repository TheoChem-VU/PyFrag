import pytest
from pyfrag.errors import PyFragSectionInputError
from pyfrag.read_input.inputblocks import InputBlocks


def test_get_input_blocks_specific_section_SCM():
    input_blocks = InputBlocks(AMS="some ams content")
    result = input_blocks.get_input_blocks_specific_section("SCM")
    assert result == {"AMS": "some ams content", "ADF": ""}


def test_get_input_blocks_specific_section_pyfrag():
    input_blocks = InputBlocks(PYFRAG="some pyfrag content")
    result = input_blocks.get_input_blocks_specific_section("PYFRAG")
    assert result == "some pyfrag content"


def test_get_input_blocks_specific_section_jobsub():
    input_blocks = InputBlocks(JOBSUB="some jobsub content")
    result = input_blocks.get_input_blocks_specific_section("JOBSUB")
    assert result == "some jobsub content"


def test_get_input_blocks_specific_section_invalid():
    input_blocks = InputBlocks()
    with pytest.raises(PyFragSectionInputError):
        input_blocks.get_input_blocks_specific_section("INVALID")  # type: ignore # it is not one of the literal options which is the expected behavior


def test_validate_missing_pyfrag():
    input_blocks = InputBlocks()
    with pytest.raises(PyFragSectionInputError, match="PyFrag section is missing in the input file"):
        input_blocks.validate()


def test_validate_missing_adf_ams():
    input_blocks = InputBlocks(PYFRAG="some pyfrag content")
    with pytest.raises(PyFragSectionInputError, match="ADF and AMS sections are missing in the input file"):
        input_blocks.validate()


def test_validate_success():
    input_blocks = InputBlocks(PYFRAG="some pyfrag content", AMS="some ams content")
    input_blocks.validate()  # Should not raise an exception


def test_iter():
    input_blocks = InputBlocks(
        JOBSUB="some jobsub content",
        PYFRAG="some pyfrag content",
        ADF="some adf content",
        AMS="some ams content",
        FRAGMENT_EXTRA={1: "fragment1 content", 2: "fragment2 content"},
        FRAGMENT_OPEN_EXTRA={1: "fragment1 open content"},
        COMPLEX_EXTRA="some complex extra content",
    )
    result = list(input_blocks)
    expected = [
        ("JOBSUB", "some jobsub content"),
        ("PYFRAG", "some pyfrag content"),
        ("ADF", "some adf content"),
        ("AMS", "some ams content"),
        ("FRAGMENT1_EXTRA", "fragment1 content"),
        ("FRAGMENT2_EXTRA", "fragment2 content"),
        ("FRAGMENT1_OPEN_EXTRA", "fragment1 open content"),
        ("COMPLEX_EXTRA", "some complex extra content"),
    ]
    assert result == expected


def test_getitem():
    input_blocks = InputBlocks(PYFRAG="some pyfrag content")
    assert input_blocks["PYFRAG"] == "some pyfrag content"


def test_str():
    input_blocks = InputBlocks(PYFRAG="some pyfrag content")
    assert "PYFRAG" in str(input_blocks)
    assert "some pyfrag content" in str(input_blocks)


def test_contains():
    input_blocks = InputBlocks(PYFRAG="some pyfrag content")
    assert "PYFRAG" in input_blocks
    assert input_blocks["JOBSUB"] == ""
