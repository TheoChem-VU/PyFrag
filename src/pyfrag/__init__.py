import logging as log

from pyfrag.environment_tools.add_scm_libbase_symlink import add_scm_libbase_to_sys_path


def initialize_pyfrag_program(log_level: int = log.INFO):
    """
    This function initializes the PyFrag program. It includes:
    - Setting up the python logger with the given log level.

    Args:
        log_level (int): The log level to set the logger to. Default is logging.INFO.
    """
    # Set up the logger again with the desired settings
    log.basicConfig(format="[%(asctime)s][%(levelname)s][%(name)s]: %(message)s", datefmt="%I:%M:%S", level=log_level)
    logger = log.getLogger(name="PyFrag Initializer")

    # Set up the enviroment for the SCM libbase, which is necessary for the correct input parsing from text to plams.Settings
    logger.info("Adding SCM libbase to sys.path...")
    add_scm_libbase_to_sys_path()
