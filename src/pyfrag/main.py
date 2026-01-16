"""
PyFrag main entry point module.

This module provides the main entry point for the PyFrag package,
handling command line arguments and dispatching to appropriate executables.
"""

import argparse
import platform
import subprocess
import sys
from pathlib import Path
from typing import Dict, Tuple, Union

from pyfrag.executables.adf.errors import ExecutableNotSupportedError, ExecutablePathNotFoundError, PyFragInputFileNotFoundError
from pyfrag.parser_factory import ExecutableType, get_parser


def print_help() -> None:
    """Print help information."""
    valid_executables = [f"'{e.value}'" for e in ExecutableType]
    help_text = f"""
PyFrag: a Python package that facilitates the analysis of reaction mechanisms in a more efficient and user-friendly way.

Usage: pyfrag <input_file> [-h] [-x <executable>] [-v]

Positional arguments:
    input_file          Path to the input file (required)

Options:
    -h, --help          Show this help message and exit
    -x, --executable    Specify the executable to run (default: 'adf')
                        Valid options: {", ".join(valid_executables)}
    -v, --verbose       Enable verbose output which determines whether the generated run script will be deleted

Examples:
    pyfrag job.in
    pyfrag job.in -x orca

See the examples folder on the GitHub page for more input file examples: https://github.com/TheoChem-VU/PyFrag.

"""
    print(help_text)


def validate_executable(value: str) -> ExecutableType:
    """Validate and convert string to ExecutableType."""
    try:
        validated_executable = ExecutableType(value.lower())
    except ValueError:
        valid_types = [e.value for e in ExecutableType]
        raise ExecutableNotSupportedError(value, valid_types)

    return validated_executable


def get_executable_path(executable_name: str) -> Union[Path, None]:
    """Get the path to a PyFrag executable."""
    executable_name = executable_name.lower()

    source_path = Path(__file__).parent.resolve()

    executable_path = source_path / "executables" / executable_name / f"{executable_name}.py"

    if not executable_path.is_file():
        raise ExecutablePathNotFoundError(f"Executable path for '{executable_name}' not found. The path was {executable_path}")

    return executable_path


def get_jobsub_section_from_input_file(input_file_sections: Dict[str, str], job_sub_key: str) -> Tuple[str, str, bool]:
    """Extracts the jobsub section from the input file sections.

    Args:
        input_file_sections (Dict[str, str]): Dictionary containing sections of the input file.

    Returns:
        Tuple[str, bool]: The jobsub section content, postambles, and a boolean indicating if it uses a job scheduler.
    """
    use_job_scheduler = False
    postamble = ""
    jobsub_content = input_file_sections.get(job_sub_key, "")

    # If there is no jobsub section, return a default one
    if not jobsub_content:
        return "!/bin/bash", "", use_job_scheduler

    # If there are certain postambles (e.g. cleanup commands), we need to account for them
    # These are recognized by the lines after the [postamble] mark which will be removed from the rest of the jobsub_content
    if "[postamble]" in jobsub_content:
        postamble = jobsub_content.split("[postamble]", 1)[1].strip()
        jobsub_content = jobsub_content.split("[postamble]", 1)[0].strip()

    # Also there is jobsub_content, but no single line starts with job schedular patterns such as #SBATCH, #PBS, etc., then it should also not be considered as using a job scheduler
    if any(line.startswith(("#SBATCH", "#PBS")) for line in jobsub_content.splitlines()):
        use_job_scheduler = True

    return jobsub_content, postamble, use_job_scheduler


def prepare_run_script(job_sub_content: str, input_file: Path, postamble: str, executable: ExecutableType) -> Path:
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

    [POSTAMBLE]
    {postamble}

    """
    job_script_file = input_file.with_suffix(".run")

    # Use .bat extension on Windows, .run elsewhere
    if platform.system() == "Windows":
        job_script_file = input_file.with_suffix(".bat")
        python_environment_line = f"python {get_executable_path(executable.value)} {input_file}"
    else:
        job_script_file = input_file.with_suffix(".run")
        if executable == ExecutableType.ADF:
            python_environment_line = f"amspython {get_executable_path(executable.value)} {input_file}"
        else:
            python_environment_line = f"python {get_executable_path(executable.value)} {input_file}"

    job_script_content = f"{job_sub_content}\n\n{python_environment_line}\n{postamble}"

    with open(job_script_file, "w") as f:
        f.write(job_script_content)

    # Ensure the job script is executable (no effect on Windows, but safe)
    try:
        job_script_file.chmod(0o755)
    except Exception:
        pass

    print(f"Run script created: {job_script_file.name}")
    return job_script_file


def parse_input() -> argparse.Namespace:
    """Main Python-based entry point for PyFrag."""
    parser = argparse.ArgumentParser(description="Python interface for the PyFrag package.", add_help=False)
    parser.add_argument("-x", "--executable", default="adf", help=f"Specify executable to run. Available options are: {', '.join([e.value for e in ExecutableType])}.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output including debug print statements.")
    parser.add_argument("input_file", nargs="?", help="Input file path")

    # Print help if no arguments are given
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)

    args, unknown_args = parser.parse_known_args()

    # Print help if input_file is missing
    if args.input_file is None:
        print_help()
        sys.exit(1)

    return args


def main():
    """Main (Python-based) entry point for PyFrag.

    It executes the following steps:

    - Parse command line arguments
    - Validate input file existence
    - Validate executable
    - Prepare and execute job script

    A separate job script is required for execution because:
        1. the correct Python environment needs to be activated or used, such as `amspython` for ADF (Amsterdam Density Functional);
        2. the job scheduler (if used) requires a specific script format and environment setup (e.g., variables denoted by #SBATCH for slurm and #PBS for PBS, and commands such as `module load`, `sbatch <script>`).
    """
    print("Starting PyFrag...")
    args = parse_input()

    # Validate input file exists
    input_file_path = Path(args.input_file).resolve()
    if not input_file_path.exists():
        raise PyFragInputFileNotFoundError(f"Input file '{input_file_path}' does not exist.")

    # Also the executable must exist, or must be supported
    validated_executable = validate_executable(args.executable)

    # Prepare the run file and execute the job script
    parser = get_parser()
    input_file_sections = parser(str(input_file_path))
    jobsub_section, postamble, use_job_scheduler = get_jobsub_section_from_input_file(input_file_sections=input_file_sections, job_sub_key="jobsub")
    run_script = prepare_run_script(jobsub_section, input_file_path, postamble, validated_executable)

    # Execute the job script if a job scheduler is used
    debug_mode = args.verbose if args.verbose else False
    if use_job_scheduler:
        if platform.system() == "Windows":
            print("Job schedulers like SLURM are not supported on Windows.")
            sys.exit(1)
        print(f"Submitting job script: {run_script.name}")
        subprocess.run(["sbatch", str(run_script)], check=True)
    else:
        print(f"Running directly: {run_script.name}")
        if platform.system() == "Windows":
            subprocess.run([str(run_script)], shell=True, check=True)
        else:
            subprocess.run(["bash", str(run_script)], check=True)
    run_script.unlink() if not debug_mode else None  # Keep the run script if using debug mode


if __name__ == "__main__":
    main()
