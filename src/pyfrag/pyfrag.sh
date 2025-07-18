#!/bin/bash

# ==================================================================
# Main shell script for PyFrag
# ------------------------------------------------------------------
# This script does the following:
# 1. It activates the virtual Python environment.
# 2. It parses the input file to the Python parser.

# Main arguments are the executable and the input file
# ==================================================================

# ==================================================================
# Helper functions
# ==================================================================

function activate_venv() {
    local venv_folder="${1:-venv}"
    source "${venv_folder}/bin/activate"
}

function print_help() {
    cat << EOF
Usage: pyfrag [-h] [-s] [-x executable] [...]

       -h            : print this information
       -x executable : start the executable named executable
                     : executable include adf, which use adf as engine (for both open and closed shell calculations)
                     : openorb, which print the real orbital energy from the open shell calculation
                     : single, which use adf as engine to do a series of SP calculation
                     : gaussian, which use gaussian as engine
                     : orca, which use orca as engine
                     : turbomole, use turbomole as engine
                     : adfold, which use ADF before 2020
                     : default executable is adf itself

The example executable is like as follow, in which job.in is job input
pyfrag job.in
or
pyfrag -x gaussian job.in
EOF
}

# ==================================================================
# Main script execution
# ==================================================================


# Print help if no arguments are provided
if [[ $# -eq 0 ]]; then
    print_help
    exit 0
fi

# Activate the virtual environment
# Get the root directory of PyFrag package
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYFRAGHOME="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Now use PYFRAGHOME for paths
VENV_PATH="$PYFRAGHOME/.venv"
if [[ -d "$VENV_PATH" ]]; then
    activate_venv "$VENV_PATH"
else
    echo "Critical: Virtual environment not found at $VENV_PATH. Please make sure that the HOSTPYFRAG environment variable is set to the root path of the PyFrag package."
    exit 1
fi

# Parse command-line arguments and call the pyfragparse.py script
while getopts ":hx:" opt; do
    case $opt in
        h)
            print_help
            exit 0
            ;;
        x)
            executable="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            print_help
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            print_help
            exit 1
            ;;
        *)
            echo "Unknown option: -$OPTARG" >&2
            print_help
            exit 1
            ;;
    esac
done

shift $((OPTIND - 1))  # Remove processed options from the argument list
# Check if the input file is provided
if [[ $# -ne 1 ]]; then
    echo "Error: Input file is required."
    print_help
    exit 1
fi

input_file_path=$(realpath "$1")
if [[ ! -f "$input_file_path" ]]; then
    echo "Error: Input file '$input_file_path' does not exist."
    exit 1
fi

# Finally, run the Python parser with the provided input file and executable
echo "Parsing $input_file_path with Python parser..."
python "main.py" "$input_file_path" -x "$executable"