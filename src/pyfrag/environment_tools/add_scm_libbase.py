import os
import sys
from pathlib import Path

EXTRA_PATHS: dict[str, Path] = {
    "SCM": Path("scripting") / "scm",
    "SCM_LIBBASE": Path("scripting") / "scm" / "libbase",
    "SCM_LIBBASE_PYTHON": Path("python3.8") / "Lib" / "site-packages"

}

def add_scm_libbase_to_sys_path():
    """Add the SCM libbase to the sys.path for Python module search."""
    AMSBIN_PATH = os.environ.get("AMSBIN")

    if AMSBIN_PATH is not None:
        AMSBIN_PATH = Path(AMSBIN_PATH)

        for description, path in EXTRA_PATHS.items():
            # Adjust the path based on the description
            full_path = AMSBIN_PATH.parent / path if description in ["SCM", "SCM_LIBBASE"] else AMSBIN_PATH / path

            if full_path.exists():
                sys.path.append(str(full_path))
            else:
                print(f"{description} path does not exist: {full_path}")

    else:
        print("AMSBIN environment variable is not set.")

add_scm_libbase_to_sys_path()

import scm  # noqa
import libbase # noqa