from pyfrag.read_input.parse_input_file import extract_section_blocks_from_file_content


def test_spellcheck_correction_fragment_extra():
    """Test the get_one_section function that returns the section name if it is supported."""
    file_content = """
    Fragment1_EXTRA
        System
            Charge 0
        End
    fraGMENt1_EXTRA END
    """
    sections_content = extract_section_blocks_from_file_content(file_content)
    assert sections_content.FRAGMENT_EXTRA == {1: "System\n            Charge 0\n        End"}


def test_spellcheck_correction_jobsub():
    """Test the get_one_section function that returns the section name if it is supported."""
    file_content = """
    job sub
    #!/bin/bash
    job sub END
    """
    sections_content = extract_section_blocks_from_file_content(file_content)
    assert sections_content.JOBSUB == "#!/bin/bash"


def test_spellcheck_correction_complex_extra():
    """Test the get_one_section function that returns the section name if it is supported."""
    file_content = """
    cOmPLEx EXTRA
    System
        Charge -1
    End
    COmPLEX EXTRA END
    """
    sections_content = extract_section_blocks_from_file_content(file_content)
    assert sections_content.COMPLEX_EXTRA == "System\n        Charge -1\n    End"


def test_extract_sections_jobsub():
    """Test the extract_sections function that extracts the content of each section from the input file."""
    file_content = """
JobSub
#!/bin/bash
#SBATCH -N 1
JOBSUB END
    """
    sections_content = extract_section_blocks_from_file_content(file_content)
    assert sections_content.JOBSUB == "#!/bin/bash\n#SBATCH -N 1"


def test_extract_sections_pyfrag():
    """Test the extract_sections function that extracts the content of each section from the input file."""
    file_content = """
PyFrag
ircpath PES1normdist57ang80.ams.amv
PyFrag END
    """
    sections_content = extract_section_blocks_from_file_content(file_content)
    assert sections_content.PYFRAG == "ircpath PES1normdist57ang80.ams.amv"


def test_extract_sections_ams():
    """Test the extract_sections function that extracts the content of each section from the input file."""
    file_content = """
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
    """
    sections_content = extract_section_blocks_from_file_content(file_content)
    assert "System" in sections_content.AMS
    assert "XC" in sections_content.AMS
    assert "BASIS" in sections_content.AMS
    assert "EndEngine" in sections_content.AMS


def test_extract_sections_fragment1_extra():
    """Test the extract_sections function that extracts the content of each section from the input file."""
    file_content = """
fragment1 EXTRA
    System
        Charge 0
    End
fragment1 EXTRA END
    """
    sections_content = extract_section_blocks_from_file_content(file_content)
    assert "System" in sections_content.FRAGMENT_EXTRA[1]
    assert "Charge 0" in sections_content.FRAGMENT_EXTRA[1]


def test_extract_sections_fragment2_extra():
    """Test the extract_sections function that extracts the content of each section from the input file."""
    file_content = """
fragment2 EXTRA
    System
        Charge -1
    End
fragment2 EXTRA END
    """
    sections_content = extract_section_blocks_from_file_content(file_content)
    assert "System" in sections_content.FRAGMENT_EXTRA[2]
    assert "Charge -1" in sections_content.FRAGMENT_EXTRA[2]


def test_extract_sections_complex_extra():
    """Test the extract_sections function that extracts the content of each section from the input file."""
    file_content = """
complex EXTRA
    System
        Charge -1
    End
complex EXTRA END
    """
    sections_content = extract_section_blocks_from_file_content(file_content)
    assert "System" in sections_content.COMPLEX_EXTRA
    assert "Charge -1" in sections_content.COMPLEX_EXTRA


def test_extract_sections_small_variations_with_comment():
    file_content_with_comment = """
PyFrag # This is a comment
    ircpath PES1normdist57ang80.ams.amv
PyFrag END
    """
    sections_content = extract_section_blocks_from_file_content(file_content_with_comment)
    assert sections_content.PYFRAG == ""


def test_extract_sections_small_variations_case_invariant():
    file_content_case_invariant = """
PyfRAG
    ircpath PES1normdist57ang80.ams.amv
PyfRAG END
    """
    sections_content = extract_section_blocks_from_file_content(file_content_case_invariant)
    assert sections_content.PYFRAG == "ircpath PES1normdist57ang80.ams.amv"


def test_extract_commented_out_jobsub_section():
    """This function should actually not pass!"""
    file_content = """
#JOBSUB
#!/bin/bash
#SBATCH -N 1

echo module list

#JOBSUB END
"""
    sections_content = extract_section_blocks_from_file_content(file_content)
    # assert sections_content.JOBSUB == "#!/bin/bash\n#SBATCH -N 1\n\necho module list\n\n#"
    assert sections_content.JOBSUB == ""
