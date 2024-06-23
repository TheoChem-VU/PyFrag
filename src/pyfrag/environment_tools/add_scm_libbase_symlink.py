import logging
import os
import sys
from pathlib import Path
from typing import Dict

import scm.plams as plams

# Define extra paths relative to AMSHOME
EXTRA_PATHS: Dict[str, Path] = {
    "SCM": Path("scripting/scm"),
    "SCM_LIBBASE": Path("scripting/scm/libbase"),
    "SCM_LIBBASE_PYTHON": Path("bin/python3.8/lib/python3.8/site-packages"),
}


def add_scm_libbase_to_sys_path():
    """Add SCM libbase to sys.path and create necessary symlinks.

    The environment variable AMSHOME must be set to the root directory of the SCM installation.
    There are two reasons for why the path needs to be updated and a symlink created:

    1. The SCM libbase is required for the correct input parsing from text to plams.Settings in the AMSJob.from_input method.
    2. If you install plams locally (via pip install), it makes a scm folder in site-packages which only contains plams and not the SCM libbase. Yet, the AMSJob.from_input method uses the SCM libbase as part in the AMS environemnt.
    So, we need to create a symlink to the SCM libbase in the LOCAL site-packages folder which links to both the scm.scripting path and the scm.libbase path (in the site-packages of SCM's python3.8 environment).

    Important folders in the SCM package:
    - Scripting folder: containing the main SCM python scripts
    - Libbase folder in the site-packages of SCM's python3.8 environment: containing the INTERNAL SCM libbase

    The result of this function is this:

    LOCAL python environment (site-packages of python x.x):
    - scm
        - plams (local copy)
        - libbase (symlink to SCM's libbase folder)

    and any imports from scm.libbase will be correctly resolved as it has also been added to sys.path so that python environment manager can find it.
    """
    # Configure logging
    AMSHOME = os.getenv("AMSHOME")
    if not AMSHOME:
        logging.error("AMSHOME environment variable is not set.")
        raise ValueError("AMSHOME environment variable is not set.")

    AMSHOME_PATH = Path(AMSHOME)
    # Append paths to sys.path
    for description, relative_path in EXTRA_PATHS.items():
        full_path = AMSHOME_PATH / relative_path
        if full_path.exists():
            sys.path.append(str(full_path))
        else:
            logging.warning(f"{description} path does not exist: {full_path}")

    # Create or update symlink for SCM libbase
    update_scm_libbase_symlink(AMSHOME_PATH)


def update_scm_libbase_symlink(amshome_path: Path):
    """Create or update the SCM libbase symlink."""
    local_scm_path = Path(plams.__file__).parent.parent  # Local SCM directory
    scm_libbase_path = amshome_path / EXTRA_PATHS["SCM_LIBBASE"]  # Target SCM libbase directory

    local_libbase_path = local_scm_path / "libbase"
    if local_libbase_path.is_symlink() or local_libbase_path.exists():
        # Remove if it's an incorrect symlink or if it somehow exists as a file/directory
        if not local_libbase_path.is_symlink() or local_libbase_path.resolve() != scm_libbase_path:
            local_libbase_path.unlink()
            local_libbase_path.symlink_to(scm_libbase_path)
            logging.info(f"Updated symlink for SCM libbase at {local_libbase_path}")
    else:
        local_libbase_path.symlink_to(scm_libbase_path)
        logging.info(f"Created symlink for SCM libbase at {local_libbase_path}")

    # Import after ensuring the symlink is correctly set up and can be used in other parts of the program
    from scm.libbase import InputParser  # noqa  # type: ignore


if __name__ == "__main__":
    add_scm_libbase_to_sys_path()
