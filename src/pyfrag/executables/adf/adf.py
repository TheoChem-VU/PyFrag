import argparse as arg
import logging
import pathlib as pl
import shutil
import sys
from pathlib import Path

from adf_driver import handle_restart, pyfrag_driver, settings_from_ams_block, write_table
from adf_to_ams_converter import main_converter
from input import process_user_input
from scm.plams import config, finish, init

sys.path.append(str(pl.Path(__file__).parent.parent.parent))  # path to the adf_to_ams_input_converter.py file

"""
Pyfrag-ADF module
Authors: Xiaobo Sun; Thomas Soini; Siebe Lekanne Deprez

This program has the following functionalities:
1: Reads in a series of Linear Transit or IRC structures (coordinate files or t21).
2: For each point PyFrag generates single point calculations for the individual fragments (based on user defined fragments)
   and the whole complex system.
3: In the following the corresponding ADF calculations are conducted.
4: The program will generate a text file containing the decomposition energies plus other, user defined, values such as the strain energy.

Example use:
amspython adf.py job_name.in
"""


def setup_logging(job_name: str, log_level: int, input_file: Path) -> logging.Logger:
    logging.basicConfig(level=log_level, filename=input_file.with_suffix(".log"), filemode="w")
    logger = logging.getLogger("PyFrag-ADF")
    logger.info(f"Starting PyFrag-ADF for job: {job_name}")
    return logger


parser = arg.ArgumentParser(description="PyFrag ADF calculations from input file")
parser.add_argument("input_file", type=Path, help="Input file containing PyFrag configuration (e.g., [job_name].in)")

args = parser.parse_args()
inputKeys = process_user_input(args.input_file)

# Set up logging
logger = setup_logging(inputKeys["job_name"], inputKeys["log_level"], args.input_file)

# Handle restart | job_name is the name of the restart directory
inputKeys["jobstate"] = handle_restart(inputKeys["job_name"])

init(folder=inputKeys["job_name"])
workdir_path = config.default_jobmanager.workdir

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

# Logging settings
for system, specific_sett in zip(["All systems", "Frag1", "Frag2", "Complex"], [settings_general, frag1_settings, frag2_settings, complex_settings]):
    logger.info(f"Settings for {system}:\n" + str(specific_sett))

# Execute PyFrag calculations: for each point on the IRC/LT, perform single point calculations for the fragments and the complex
tableValue, inputKeys = pyfrag_driver(inputKeys, frag1_settings, frag2_settings, complex_settings)

logger.info("Writing table to file and removing extra files")
write_table(tableValue, inputKeys["job_name"])

# Remove the restart directory if it exists.
if inputKeys["jobstate"] is not None:
    shutil.rmtree(inputKeys["jobstate"])
finish()
logging.info("PyFrag finished")
