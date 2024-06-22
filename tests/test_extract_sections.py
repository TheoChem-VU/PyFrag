import pytest
from pyfrag.errors import UnsupportedSectionError
from pyfrag.input.parse_input_file import extract_section_blocks_from_file_content, get_one_section


def test_get_one_section():
    """Test the get_one_section function that returns the section name if it is supported."""
    assert get_one_section("JOBSUB") == "JOBSUB"
    assert get_one_section("NOT A SECTION") is None
    assert get_one_section("PYFRAG") == "PYFRAG"
    with pytest.raises(UnsupportedSectionError):
        get_one_section("JOB SUB")
    with pytest.raises(UnsupportedSectionError):
        get_one_section("FRAGMENT1_EXTRA")
    with pytest.raises(UnsupportedSectionError):
        get_one_section("COMPLEX_EXTRA")


def test_extract_sections():
    """Test the extract_sections function that extracts the content of each section from the input file."""
    file_content = """
JobSub
#!/bin/bash
#SBATCH -N 1
JOBSUB END

PyFrag
ircpath PES1normdist57ang80.ams.amv
PyFrag END

AMS
    System
        Charge 0
    End
    Task SinglePoint
    Engine ADF
        XC
            GGA OLYP
        END
        NumericalQuality normal
        BASIS
            TYPE DZP
            CORE None
        END
    EndEngine
AMS END

fragment1 EXTRA
    System
        Charge 0
    End
fragment1 EXTRA END

fragment2 EXTRA
    System
        Charge -1
    End
fragment2 EXTRA END

complex EXTRA
    System
        Charge -1
    End
complex EXTRA END
    """
    sections_content = extract_section_blocks_from_file_content(file_content)

    assert sections_content["JOBSUB"] == "#!/bin/bash\n#SBATCH -N 1"
    assert sections_content["PYFRAG"] == "ircpath PES1normdist57ang80.ams.amv"

    # Asserting word for word instead of the whole string because of the indentation
    assert "System" in sections_content["AMS"]
    assert "XC" in sections_content["AMS"]
    assert "BASIS" in sections_content["AMS"]
    assert "EndEngine" in sections_content["AMS"]

    assert "System" in sections_content["FRAGMENT1 EXTRA"]
    assert "Charge 0" in sections_content["FRAGMENT1 EXTRA"]

    assert "System" in sections_content["FRAGMENT2 EXTRA"]
    assert "Charge -1" in sections_content["FRAGMENT2 EXTRA"]

    assert "System" in sections_content["COMPLEX EXTRA"]
    assert "Charge -1" in sections_content["COMPLEX EXTRA"]
