#!/bin/bash
# PyFrag Virtual Environment Activation Script
# This script ensures the virtual environment is activated and the correct Python is used

# Function to find and activate the PyFrag virtual environment
activate_pyfrag_venv() {
    local pyfrag_home="${HOSTPYFRAG:-${PYFRAGHOME}}"

    if [ -z "$pyfrag_home" ]; then
        echo "Error: HOSTPYFRAG or PYFRAGHOME environment variable must be set" >&2
        exit 1
    fi

    # Check for virtual environment
    local venv_path="$pyfrag_home/.venv"

    if [ -d "$venv_path" ] && [ -f "$venv_path/bin/activate" ]; then
        echo "Activating PyFrag virtual environment: $venv_path"
        source "$venv_path/bin/activate"

        # Verify activation worked
        if [ -n "$VIRTUAL_ENV" ]; then
            echo "Virtual environment activated successfully"
            echo "Python executable: $(which python)"
            echo "Python version: $(python --version)"
        else
            echo "Warning: Virtual environment activation may have failed"
        fi
    else
        echo "Warning: Virtual environment not found at $venv_path"
        echo "Using system Python: $(which python3 || which python)"
    fi
}

# Function to run PyFrag with proper environment
run_pyfrag_cluster() {
    activate_pyfrag_venv

    # Use python from virtual environment if available, otherwise fall back to python3/python
    local python_cmd
    if [ -n "$VIRTUAL_ENV" ]; then
        python_cmd="python"
    else
        python_cmd=$(which python3 || which python || echo "python")
    fi

    echo "Executing: $python_cmd $@"
    exec "$python_cmd" "$@"
}

# Export the function so it can be used in sub scripts
export -f activate_pyfrag_venv
export -f run_pyfrag_cluster
