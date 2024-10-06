from pyfrag.enums import CalcType
from pyfrag.read_input.read_ams_section import determine_calc_type


def test_determine_calc_type_unrestricted_no_spinpolarization():
    """Test the get_one_section function that returns the section name if it is supported."""
    file_content = """
    System
    Charge 0
    End

    Engine ADF
    FragOccupations
        CX
        A1 8 // 7
        A2 0 // 0
        B1 2 // 2
        B2 4 // 3
        subend
        NH2
        A1 3 // 4
        A2 1 // 1
        B1 1 // 1
        B2 3 // 4
        subend
    End

    IrrepOccupations
        A1 18 1
        A2 0
        B1 6
        B2 8 1
    End
    EndEngine


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
    """
    calctype = determine_calc_type(file_content)
    assert calctype == CalcType.UNRESTRICTED_NO_SPINPOL


def test_determine_calc_type_unrestricted_spinpolarization():
    """Test the get_one_section function that returns the section name if it is supported."""
    file_content = """
    System
        Charge 0
    End

    Engine ADF
        UnrestricedFragments
        Unrestricted

        Spinpolarization 1


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
    """
    assert determine_calc_type(file_content) == CalcType.UNRESTRICTED_SPINPOL


def test_determine_calc_type_restricted():
    """Test the get_one_section function that returns the section name if it is supported."""
    file_content = """
    System
    Charge 0
    End

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
    """
    assert determine_calc_type(file_content) == CalcType.RESTRICTED
