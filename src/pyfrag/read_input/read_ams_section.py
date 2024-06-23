import logging

from pyfrag.process_input.enums import CalcType

logger = logging.getLogger(name="AMS/ADF Section Reader")


def determine_calc_type(ams_section: str) -> CalcType:
    """
    Determines the type of calculation from the ADF input file section. This includes the optional Fragment en Complex blocks.
    It searches for the unrestricted keywords, and for the fragoccupations keywords.
    """
    ams_section = ams_section.lower()

    if any(x in ams_section for x in ["fragoccupations", "irrepoccupations"]):
        return CalcType.UNRESTRICTED_NO_SPINPOL

    if any(x in ams_section for x in ["unrestricted", "spinpolarization"]):
        return CalcType.UNRESTRICTED_SPINPOL

    return CalcType.RESTRICTED
