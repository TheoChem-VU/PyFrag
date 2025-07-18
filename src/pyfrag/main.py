"""
PyFrag main entry point module.

This module provides the main entry point for the PyFrag package,
handling command line arguments and dispatching to appropriate executables.
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Dict, Tuple, Union

from pyfrag.config import get_config
from pyfrag.parser_factory import ExecutableType, get_parser_from_executable


def get_executable_path(executable_name: str) -> Union[Path, None]:
    """Get the path to a PyFrag executable."""
    executable_name = executable_name.lower()
    config = get_config()

    host_bin_path = config.get_executables_path(executable_name)
    if host_bin_path.exists():
        return host_bin_path

    return None


def validate_executable(value: str) -> ExecutableType:
    """Validate and convert string to ExecutableType."""
    try:
        return ExecutableType(value.lower())
    except ValueError:
        valid_types = [e.value for e in ExecutableType]
        raise argparse.ArgumentTypeError(f"Invalid executable '{value}'. Valid options are: {', '.join(valid_types)}")


def get_jobsub_section_from_input_file(input_file_sections: Dict[str, str]) -> Tuple[str, bool]:
    """Extracts the jobsub section from the input file sections.

    Args:
        input_file_sections (Dict[str, str]): Dictionary containing sections of the input file.

    Returns:
        Tuple[str, bool]: The jobsub section content and a boolean indicating if it uses a job scheduler.
    """
    jobsub_key = "jobsub"
    jobsub_content = input_file_sections.get(jobsub_key, "")

    if not jobsub_content:
        return "!/bin/bash", False

    return jobsub_content, True


def prepare_run_script(job_sub_content: str, input_file: Path, venv_path: Union[Path, None], executable: ExecutableType) -> Path:
    """Makes a .run file for the job script and runs it.

    Example job script content:
    [JOBSUB SECTION]
    If use SLURM or similar job scheduler, the job script might look like this:
        #!/bin/bash
        #SBATCH --job-name=pyfrag_job
        #SBATCH --output=pyfrag_job.out
        #SBATCH --error=pyfrag_job.err
        #SBATCH --time=01:00:00
        #SBATCH --ntasks=1
        #SBATCH --cpus-per-task=4

    else no job scheduler is used:
        #!/bin/bash

    [PYTHON ENVIRONMENT SECTION]
    if ADF, use amspython (for the InputParser that calls the compiled libbase InputParser of the AMS package):
        amspython /path/to/pyfrag/[executable].py input_file.in
    else use the virtual Python environment:
        source /path/to/venv/bin/activate  # Activate virtual environment if needed
        python /path/to/pyfrag/[executable].py input_file.in

    """
    job_script_file = input_file.with_suffix(".run")

    if executable == ExecutableType.ADF:
        python_environment_line = f"amspython {get_executable_path(executable.value)} {input_file}"
    else:
        if venv_path:
            python_environment_line = f"source {venv_path}/bin/activate && python {get_executable_path(executable.value)} {input_file}"
        else:
            python_environment_line = f"python {get_executable_path(executable.value)} {input_file}"

    job_script_content = f"{job_sub_content}\n\n{python_environment_line}\n"

    with open(job_script_file, "w") as f:
        f.write(job_script_content)

    # Ensure the job script is executable
    job_script_file.chmod(0o755)

    print(f"Run script created: {job_script_file}")
    return job_script_file


def print_help() -> None:
    """Print help information."""
    valid_executables = ", ".join([e.value for e in ExecutableType])
    help_text = f"""
Python interface for the PyFrag package.

Usage: pyfrag <input_file> [options]

The script automatically detects cluster mode if a JOBSUB section is present in the input file.

Positional arguments:
    input_file          Path to the input file (required)

Options:
    -h, --help          Show this help message and exit
    -x, --executable    Specify the executable to run (default: 'adf')
                        Valid options: {valid_executables}

Examples:
    # Direct mode (no JOBSUB section)
    pyfrag job.in
    pyfrag job.in -x gaussian

    # Cluster mode (JOBSUB section present in input file)
    pyfrag job_with_jobsub.in
    pyfrag job_with_jobsub.in -x adf

"""
    print(help_text)


def main():
    """Main Python-based entry point for PyFrag."""
    parser = argparse.ArgumentParser(add_help=False)  # We'll handle help ourselves
    parser.add_argument("-h", "--help", action="store_true", help="Show help message")
    parser.add_argument("-x", "--executable", default="adf", help=f"Specify executable to run. Available options are {list(ExecutableType)}.")
    parser.add_argument("input_file", help="Input file path")

    args, unknown_args = parser.parse_known_args()
    validated_executable = validate_executable(args.executable)

    if args.help:
        print_help()
        sys.exit(0)

    # Validate input file exists
    input_file_path = Path(args.input_file)
    if not input_file_path.exists():
        print(f"Error: Input file '{args.input_file}' does not exist", file=sys.stderr)
        sys.exit(1)

    # Prepare the run file and execute the job script
    parser = get_parser_from_executable(validated_executable)
    input_file_sections = parser(str(input_file_path))
    jobsub_section, use_job_scheduler = get_jobsub_section_from_input_file(input_file_sections)
    venv_path = get_config().virtual_env

    prepare_run_script(jobsub_section, input_file_path, venv_path, validated_executable)

    # Execute the job script if a job scheduler is used
    if use_job_scheduler:
        print(f"Submitting job script: {input_file_path.with_suffix('.run')}")
        subprocess.run(["sbatch", str(input_file_path.with_suffix(".run"))], check=True)
    else:
        print(f"Running directly: {input_file_path.with_suffix('.run')}")
        subprocess.run([str(input_file_path.with_suffix(".run"))], check=True)


if __name__ == "__main__":
    main()
