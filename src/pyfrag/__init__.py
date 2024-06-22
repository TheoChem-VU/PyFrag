import logging as log
from typing import Literal


def initialize_pyfrag_program(log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO"):
    """
    This function initializes the PyFrag program. It includes:
    - Setting up the python logger with the given log level.

    Args:
        log_level (str): The log level that will be used by the logger. Options are "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL".
    """
    # Set up the logger
    log.basicConfig(format="[%(asctime)s][%(levelname)s][%(module)s]: %(message)s", datefmt="%I:%M:%S", level=log.WARNING)
