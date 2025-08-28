import argparse as arg
import logging
import shutil
import sys
from pathlib import Path
from typing import Tuple

from adf_driver import handle_restart, pyfrag_driver, settings_from_ams_block, write_table
from adf_to_ams_converter import main_converter
from input import InputKeys, process_user_input
from scm.plams import Settings, config, finish, init

"""
Pyfrag-ADF module
Authors: Xiaobo Sun; Thomas Soini; Siebe Lekanne Deprez

This program has the following functionalities:
1: Reads in a coordinate file containing, for example, an Intrinsic Reaction Coordinate trajectory. Supported are xyz, ams.rkf, and amv formats.
2: For each point on the trajectory, PyFrag performs a fragment analysis calculation (single point calculations) for the individual fragments (based on user defined fragments)
   and the whole complex system.
3: The program extracts relevant data that the user defines from the ADF calculations.
4: The program generates a text file containing the decomposition energies plus other, user defined, values such as orbital overlaps.

Example use:
amspython adf.py job_name.in
"""


def pyfrag_program_string() -> str:
    """Returns a string that is printed at the start of the program"""
    return """
===============================================================================
                              PyFrag - ADF
===============================================================================
Program for performing fragment analysis calculations to facilitate Activation-
Strain Analyses (ASM/ASE) and Energy Decomposition Analyses (EDA).

For any questions or feedback, please contact the developers via Github:
https://github.com/TheoChem-VU/PyFrag
===============================================================================
"""


def setup_logging(job_name: str, log_level: int) -> logging.Logger:
    logging.basicConfig(level=log_level, filemode="w", stream=sys.stdout)
    logger = logging.getLogger("PyFrag-ADF")
    logger.info(f"Starting PyFrag-ADF for job: {job_name}")
    return logger


def extract_plams_settings_from_adf_settings(inputKeys: InputKeys) -> Tuple[Settings, ...]:
    # This (ugly) block of code is necessary to check if the user has provided an AMS input file or an old ADF input file
    old_ADF_input = False
    if inputKeys["adfinputfile"] is not None:
        settings_general = settings_from_ams_block(inputKeys["adfinputfile"])
    elif inputKeys["old_adfinputfile"] is not None:
        old_ADF_input = True
        settings_general = main_converter(inputKeys["old_adfinputfile"])
    else:
        raise ValueError("Failed to find the AMS input options in the input file.\nPlease make sure that you have the 'AMS / AMS End' or 'ADF / ADF END' blocks in the PyFrag input file.")

    frag1_settings = settings_general.copy()  # copy is necessary to avoid linking errors when changing settings
    frag2_settings = settings_general.copy()
    complex_settings = settings_general.copy()

    # Update (and possibly convert the <2019 ADF parsing to >2019 AMS) settings for the fragments and the complex
    for extra_input, extra_settings in zip(["fragment1_extra", "fragment2_extra", "complex_extra"], [frag1_settings, frag2_settings, complex_settings]):
        if inputKeys[extra_input] is None:
            continue

        extra_input = inputKeys[extra_input]

        if old_ADF_input:
            extra_settings.update(main_converter(extra_input))
        else:
            extra_settings.update(settings_from_ams_block(extra_input))

    return frag1_settings, frag2_settings, complex_settings


def main():
    # =====================================================================
    # Parse the inputs and set up logging for the ADF calculations
    # =====================================================================

    print(pyfrag_program_string())
    parser = arg.ArgumentParser(description="PyFrag ADF calculations from input file")
    parser.add_argument("input_file", type=Path, help="Input file containing PyFrag configuration (e.g., [job_name].in)")

    args = parser.parse_args()
    inputKeys = process_user_input(args.input_file)

    # Set up logging and print the extracted input from the PyFrag input file (parsed by the user)
    logger = setup_logging(inputKeys["job_name"], inputKeys["log_level"])
    logger.info("Processed input keys:")
    for key, value in inputKeys.items():
        if not value:
            continue  # Skip empty values

        if key in ["adfinputfile", "old_adfinputfile", "fragment1_extra", "fragment2_extra", "complex_extra"]:
            logger.info(f" - {key}")
        else:
            logger.info(f" - {key}: {value}")

    # =====================================================================
    # Handle restart directory and initialization
    # =====================================================================

    inputKeys["restart_dir_name"] = handle_restart(inputKeys["job_name"])

    # If the folder that will be created by plams already exists, remove it first in order to prevent [job_name].xxx
    if (args.input_file.parent / inputKeys["job_name"]).is_dir():
        shutil.rmtree(args.input_file.parent / inputKeys["job_name"])

    init(folder=inputKeys["job_name"])
    config.log.file = 5  # Quite verbose logging in the log file (7 is most verbose) created in the plams folder. This is for better debugging purposes
    config.log.stdout = 0  # No plams logs to the stdout because only PyFrag-related logs should be shown

    frag1_settings, frag2_settings, complex_settings = extract_plams_settings_from_adf_settings(inputKeys)

    # Logging settings
    for system, specific_sett in zip(["Frag1", "Frag2", "Complex"], [frag1_settings, frag2_settings, complex_settings]):
        logger.info(f"Settings for {system}:\n" + str(specific_sett))

    # =====================================================================
    # Execute the PyFrag calculations and write the results to a table:
    # perform single point calculations for each point on the IRC/LT for the fragments and the complex
    # =====================================================================

    tableValue, inputKeys = pyfrag_driver(inputKeys, frag1_settings, frag2_settings, complex_settings)

    logger.info("Writing table to file and removing extra files")
    write_table(tableValue, inputKeys["job_name"])

    # =====================================================================
    # Clean up and finish
    # =====================================================================

    if inputKeys["restart_dir_name"] is not None:
        shutil.rmtree(inputKeys["restart_dir_name"])
    finish()

    logging.info("PyFrag finished")


if __name__ == "__main__":
    main()
