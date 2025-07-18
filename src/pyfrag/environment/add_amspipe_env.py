import os
import shutil
import site
import sys
from pathlib import Path


def find_scm_site_packages():
    for path in site.getsitepackages() + [site.getusersitepackages()]:
        scm_path = Path(path) / "scm"
        if scm_path.is_dir():
            return scm_path
    raise FileNotFoundError("Could not find 'scm' directory in site-packages.")


def copy_folder_to_scm(src_folder):
    scm_path = find_scm_site_packages()
    print(f"Found 'scm' directory at: {scm_path}")
    src_folder = Path(src_folder)
    dest_folder = scm_path / src_folder.name
    if dest_folder.exists():
        print(f"Destination folder '{dest_folder}' already exists. Overwriting.")
        shutil.rmtree(dest_folder)
    shutil.copytree(src_folder, dest_folder)
    print(f"Copied '{src_folder}' to '{dest_folder}'.")


if __name__ == "__main__":
    ams_home_path = os.environ.get("AMSHOME")

    if not ams_home_path:
        print("Environment variable 'AMSHOME' is not set.")
        sys.exit(1)

    ams_pipe_folders = Path(ams_home_path) / "scripting" / "scm" / "amspipe"
    if not ams_pipe_folders.is_dir():
        print(f"Source folder '{ams_pipe_folders}' does not exist or is not a directory.")
        sys.exit(1)
    copy_folder_to_scm(ams_pipe_folders)
