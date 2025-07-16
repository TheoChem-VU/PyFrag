#!/usr/bin/env python3
"""
PyFrag Cluster Entry Point

This script serves as the main entry point for PyFrag calculations on the cluster.
It handles virtual environment activation and dispatches to the appropriate
calculation backend (adf_new, gaussian, orca, etc.).
"""

import argparse
import os
import sys
from pathlib import Path


def setup_pyfrag_environment():
    """Set up the PyFrag environment and virtual environment."""
    # Get the PyFrag home directory
    script_dir = Path(__file__).parent.absolute()
    pyfrag_home = script_dir.parent.parent  # Go up from src/pyfrag/ to PyFrag root

    # Set environment variables
    os.environ["PYFRAGHOME"] = str(pyfrag_home)
    os.environ["HOSTPYFRAG"] = str(pyfrag_home)

    # Add src to Python path
    src_dir = pyfrag_home / "src"
    if str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))

    # Try to activate virtual environment if it exists
    venv_path = pyfrag_home / ".venv"
    if venv_path.exists():
        # Add virtual environment site-packages to path
        venv_site_packages = venv_path / "lib" / f"python{sys.version_info.major}.{sys.version_info.minor}" / "site-packages"
        if venv_site_packages.exists() and str(venv_site_packages) not in sys.path:
            sys.path.insert(0, str(venv_site_packages))

        # Update environment variables
        os.environ["VIRTUAL_ENV"] = str(venv_path)
        venv_bin = str(venv_path / "bin")
        current_path = os.environ.get("PATH", "")
        if venv_bin not in current_path:
            os.environ["PATH"] = f"{venv_bin}:{current_path}"


def run_backend(backend, job_dir, input_file, extra_args):
    """Run the specified PyFrag backend."""
    setup_pyfrag_environment()

    # Import PyFrag config after environment setup
    try:
        from pyfrag.config import get_config

        config = get_config()
    except ImportError as e:
        print(f"Error: Could not import PyFrag package: {e}", file=sys.stderr)
        print("Make sure PyFrag is properly installed in the virtual environment.", file=sys.stderr)
        sys.exit(1)

    # Map backend names to their standalone directories
    backend_map = {"adf": "adf_new", "adf_new": "adf_new", "adf_old": "adf_old", "adfold": "adf_old", "gaussian": "gaussian", "orca": "orca", "turbomole": "turbomole", "open": "adf_open", "adf_open": "adf_open"}

    if backend not in backend_map:
        print(f"Error: Unknown backend '{backend}'", file=sys.stderr)
        print(f"Available backends: {', '.join(backend_map.keys())}", file=sys.stderr)
        sys.exit(1)

    standalone_dir = backend_map[backend]
    pyfrag_script = config.get_host_path() / "standalone" / standalone_dir / "PyFrag.py"

    if not pyfrag_script.exists():
        print(f"Error: PyFrag script not found: {pyfrag_script}", file=sys.stderr)
        sys.exit(1)

    # Change to job directory
    original_dir = os.getcwd()
    os.chdir(job_dir)

    try:
        # Import and run the PyFrag script as a module
        print(f"Running PyFrag backend: {backend} ({standalone_dir})")
        print(f"Job directory: {job_dir}")
        print(f"Input file: {input_file}")

        # Add the standalone directory to Python path
        standalone_path = str(pyfrag_script.parent)
        if standalone_path not in sys.path:
            sys.path.insert(0, standalone_path)

        # Import the PyFrag module from the standalone directory
        import importlib.util

        spec = importlib.util.spec_from_file_location("pyfrag_backend", pyfrag_script)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not load PyFrag script from {pyfrag_script}")

        pyfrag_module = importlib.util.module_from_spec(spec)

        # Set up sys.argv as if we called the script directly
        original_argv = sys.argv.copy()
        sys.argv = ["PyFrag.py"] + extra_args

        try:
            # Execute the PyFrag script
            spec.loader.exec_module(pyfrag_module)
        finally:
            # Restore original argv
            sys.argv = original_argv

    except Exception as e:
        print(f"Error running PyFrag backend: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)
    finally:
        # Restore original directory
        os.chdir(original_dir)


def main():
    parser = argparse.ArgumentParser(description="PyFrag Cluster Entry Point")
    parser.add_argument("--backend", "-b", required=True, help="PyFrag backend to use (adf, gaussian, orca, etc.)")
    parser.add_argument("--job-dir", "-d", required=True, help="Job directory where calculation should run")
    parser.add_argument("--input-file", "-i", required=True, help="Input file for the calculation")

    # Parse known args to allow passing extra arguments to the backend
    args, extra_args = parser.parse_known_args()

    # Add the input file to extra args if not already there
    if args.input_file not in extra_args:
        extra_args.append(args.input_file)

    run_backend(args.backend, args.job_dir, args.input_file, extra_args)


if __name__ == "__main__":
    main()
