from typing import Tuple

from scm.plams import AMSJob, Molecule, Settings

from pyfrag.enums import CalcType
from pyfrag.errors import PyFragSectionInputError
from pyfrag.read_input.inputblocks import InputBlocks


def convert_input_block_to_settings(input_block: str) -> Settings:
    """
    Returns a Settings object from an input block that specifies the calculation settings for the AMS program & ADF engine.

    Adapted from the AMSJob.from_inputfile method to read in the settings from a file and return a Settings object
    Reason for adapting was because the AMSJob.from_inputfile method does remove the ams block from the settings object
    """
    try:
        pre_calc_job: AMSJob = AMSJob.from_input(input_block)
    except Exception as e:
        raise PyFragSectionInputError(f"Error occurred when converting input block to settings: {e}")

    # Validate the settings: if there is an system block (containing symmetrize or symmetry), then it is removed since it crashes the program (e.g., there is no information about "symmetry" in the settings so we can't extract the appropriate point group)
    if pre_calc_job.settings.input.ams.get("System") is not None:
        _ = pre_calc_job.settings.input.ams.pop("System") if "System" in pre_calc_job.settings.input.ams else []

    settings = pre_calc_job.settings

    # Add "Task" block if it is not present
    if "Task" not in settings.input.ams:
        settings.input.ams.Task = "SinglePoint"

    # Charge is not part of the settings, but absorbed into the molecule properties. So we need to put it back into the plams settings object.
    if pre_calc_job.molecule is not None:
        molecule: Molecule = pre_calc_job.molecule[""]  # type: ignore
        charge = molecule.properties.charge
        settings.input.ams.System.Charge = charge

    return settings


def determine_calculation_type_from_ams_input(ams_input: InputBlocks) -> CalcType:
    """
    Determine the type of calculation that is specified in the input blocks.
    Returns an Enum value that corresponds to the type of calculation.
    """
    input_blocks_content = [input_block for _, input_block in ams_input.get_input_blocks_specific_section("SCM").items()]
    input_blocks_content = [block.lower() for block in input_blocks_content if block is not None]

    if any("spinpolarization" in block for block in input_blocks_content):
        if any("fragoccupations" in block for block in input_blocks_content):
            return CalcType.UNRESTRICTED_NO_SPINPOL
        else:
            return CalcType.UNRESTRICTED_SPINPOL

    return CalcType.RESTRICTED


def process_ams_input(ams_input: InputBlocks) -> Tuple[Settings, CalcType]:
    """
    Process the input blocks that specifies settings for the AMS/ADF program.
    Returns a Settings object for each input block that can be used to run an AMSJob.
    Also returns the information regarding the type of calculation (Enum) for branching the workflow.
    """

    calc_type = determine_calculation_type_from_ams_input(ams_input)

    key = "AMS" if "AMS" in ams_input else "ADF"
    if key in ams_input and ams_input[key] is not None:
        processed_settings = convert_input_block_to_settings(ams_input[key])

    return processed_settings, calc_type
