from typing import Dict, Tuple

from scm.plams import AMSJob, Molecule, Settings

from pyfrag.config.config import pyfrag_config
from pyfrag.process_input.enums import CalcType
from pyfrag.read_input.inputblocks import InputBlocks


def convert_inputblock_to_settings(input_block: str) -> Settings:
    """
    Returns a Settings object from an inputfile

    Adapted from the AMSJob.from_inputfile method to read in the settings from a file and return a Settings object
    Reason for adapting was because the AMSJob.from_inputfile method does remove the ams block from the settings object
    """
    pre_calc_job: AMSJob = AMSJob.from_input(input_block)

    # This does not include the "System" block; it does include the "Task" and "Engine" block
    settings = pre_calc_job.settings

    # Add "Task" block if it is not present
    if "Task" not in settings.input.ams:
        settings.input.ams.Task = "SinglePoint"

    # # First make sure that the inputfile is read in correctly (parsing the heredoc)
    # with open(inputfile, 'r') as f:
    #     inp_file = parse_heredoc(f.read(), 'eor')

    # # Then read in the "System" block from the inputfile and add it to the settings object
    # unnested_settings = InputParser().to_settings("ams", inp_file)
    # settings.input["ams"].update(unnested_settings["ams"])

    # Work around for the fact that the ams block is not included in the settings object
    # The above code apparently does not allow for overriding the molecule
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
    input_blocks_content = [input_block for _, input_block in ams_input.get_input_blocks_specific_section("AMS").items()]
    input_blocks_content = [block.lower() for block in input_blocks_content if block is not None]

    if any("spinpolarization" in block for block in input_blocks_content):
        if any("fragoccupations" in block for block in input_blocks_content):
            return CalcType.UNRESTRICTED_NO_SPINPOL
        else:
            return CalcType.UNRESTRICTED_SPINPOL

    return CalcType.RESTRICTED


def process_ams_input(ams_input: InputBlocks) -> Tuple[Dict[str, Settings], CalcType]:
    """
    INTERFACE FUNCTION

    Process the input blocks that specifies settings for the AMS/ADF program.
    Returns a Settings object for each input block that can be used to run an AMSJob.
    Also returns the information regarding the type of calculation (Enum) for branching the workflow.
    """
    processed_settings: Dict[str, Settings] = {}

    calc_type = determine_calculation_type_from_ams_input(ams_input)
    ams_keys = pyfrag_config.ams.input_blocks

    for key in ams_keys:
        if key in ams_input:
            if ams_input[key] is not None:
                processed_settings[key] = convert_inputblock_to_settings(ams_input[key])

    return processed_settings, calc_type
