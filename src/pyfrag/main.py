#!/usr/bin/env python3
"""
PyFrag main entry point module.

This module provides the main entry point for the PyFrag package,
handling command line arguments and dispatching to appropriate executables.
"""

import argparse
import subprocess
import sys

from pyfrag.config import get_config


def get_executable_path(executable_name):
    """Get the path to a PyFrag executable."""
    config = get_config()

    # First check in host/bin (new location)
    host_bin_path = config.get_executables_path() / f"{executable_name}.sh"
    if host_bin_path.exists():
        return host_bin_path

    # Then check in main bin directory
    main_bin_path = config.get_scripts_path() / f"{executable_name}.sh"
    if main_bin_path.exists():
        return main_bin_path

    return None


def print_help():
    """Print help information."""
    help_text = """
Usage: pyfrag [-h] [-s] [-x command] [...]

       -h          : print this information
       -x command  : start the executable named command
                   : command include adf, which use adf as engine (for both open and closed shell calculations)
                   : openorb, which print the real orbital energy from the open shell calculation
                   : single, which use adf as engine to do a series of SP calculation
                   : gaussian, which use gaussian as engine
                   : orca, which use orca as engine
                   : turbomole, use turbomole as engine
                   : adfold, which use ADF before 2020
                   : default command is adf itself

The example command is like as follow, in which job.in is job input
pyfrag job.in
or
pyfrag -x gaussian job.in
"""
    print(help_text)


def main():
    """Main entry point for PyFrag."""
    parser = argparse.ArgumentParser(add_help=False)  # We'll handle help ourselves
    parser.add_argument("-h", "--help", action="store_true", help="Show help message")
    parser.add_argument("-x", "--executable", default="adf", help="Specify executable to run")
    parser.add_argument("-s", "--silence", action="store_true", help="Run silently")
    parser.add_argument("args", nargs="*", help="Arguments to pass to the executable")

    args, unknown_args = parser.parse_known_args()

    if args.help:
        print_help()
        sys.exit(0)

    # Set up PyFrag environment
    config = get_config()
    env = config.setup_environment()

    # Get the executable path
    executable_path = get_executable_path(args.executable)

    if not executable_path:
        # Fallback to pyfrag.sh
        fallback_path = config.get_scripts_path() / "pyfrag.sh"
        if fallback_path.exists():
            executable_path = fallback_path
        else:
            print(f"Error: Cannot find executable for '{args.executable}'", file=sys.stderr)
            sys.exit(1)

    # Prepare command
    cmd = [str(executable_path)] + args.args + unknown_args

    print(f"Running command: {' '.join(cmd)}", file=sys.stderr)

    try:
        if args.silence:
            # Run silently
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True, env=env)
        else:
            # Run normally
            subprocess.run(cmd, check=True, env=env)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)
    except FileNotFoundError:
        print(f"Error: Cannot execute {executable_path}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
