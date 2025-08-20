#!/bin/bash

# ==================================================================
# This script is used to deploy the PyFrag application on a High performance computer (HPC) cluster.
# The following is required for deploying the PyFrag application:
#
#   1. Install the Python package manager "uv" if it is not already installed
#   2. Create a virtual environment (Python) with the pyfrag package installed through `uv sync`
#   3. Verify the PyFrag installation + cleaning up
# ==================================================================

echo "Trying to deploy PyFrag on the cluster which involves:"
echo "1. Installing the Python package manager 'uv'"
echo "2. Creating a virtual environment (Python) with the pyfrag package installed through 'uv sync'"
echo "3. Uninstalling uv"


# Determine pyfrag ROOTFOLDER:
PYFRAGROOT="$(pwd)"

# NOTE: this function assumes the Linux/MacOS operating system
# We do not include yet support for Windows.
install_uv() {
    # Try first with the curl command
    curl -LsSf https://astral.sh/uv/install.sh | sh

    # Otherwise, try with wget
    wget -qO- https://astral.sh/uv/install.sh | sh

    echo "uv installed successfully"
}

# Check if activating the virtual environment is successful which should be located in the $PYFRAGROOT/.venv/bin folder
# The PyFrag package should be installed in the lib/site-packages folder, which can be checked by running the "import pyfrag" command
function verify_pyfrag_installation() {
    echo "Verifying PyFrag installation..."

    # Try to load the virtual environment
    if [ -f ".venv/bin/activate" ]; then
        source .venv/bin/activate
    else
        echo "Critical: Virtual environment not found. Installation with uv went wrong. Did you run `uv sync --no-cache --no-dev --no-editable` in the root folder of the PyFrag package?"
        exit 1
    fi

    # Check if python is available
    if ! command -v python &> /dev/null; then
        echo "Critical: Python is not available in your environment, which means the virtual environment was not activated successfully."
        exit 1
    fi

    # Check if pyfrag is importable
    python -c "import pyfrag" 2>/dev/null
    if [[ $? -ne 0 ]]; then
        echo "Critical: The 'pyfrag' package is not importable in this Python environment."
        echo "Information about the current Python environment:"
        echo "Folder: $(which python)"
        echo "Version: $(python -c "import sys; print(sys.version)")"

        echo -e "Please ensure you have installed pyfrag in your current Python environment. This can be achieved by:\n1. Changing to the PyFrag root folder\n2. Activating the Python environment (\`source .venv/bin/activate\`)\n3. Running 'pip install .'."
        exit 1
    fi
}


# 1. Install the Python package manager "uv"
if ! command -v uv &> /dev/null
then
    echo "Installing uv..."
    install_uv
else
    echo "uv is already installed."
fi

# 2. Create a virtual environment (Python) with the pyfrag package installed through `uv sync`
# Do not use cache (as this may give unexpected errors when solving dependencies)
# Do not use developer dependencies as we don't want to include dependencies for testing (pytest) or documentation (sphinx)
echo "Creating virtual environment with dependencies and installing the python package..."

if [ -d ".venv" ]; then
    rm -rf .venv
fi

uv sync --no-cache --no-dev --no-editable

# 3. Check if the virtual environment is created and can be activated
if [ -d ".venv" ]; then
    echo "Virtual environment created successfully."
else
    echo "Failed to create virtual environment."
    exit 1
fi

# 3. Uninstalling uv
# See https://docs.astral.sh/uv/getting-started/installation/#uninstallation (accessed August 18 2025)
echo "Uninstalling uv..."
uv cache clean
rm ~/.local/bin/uv ~/.local/bin/uvx

# Verify PyFrag installation
verify_pyfrag_installation

echo "Finally, please make sure to set the PYFRAGROOT environment variable in the module file to the bin folder so that the 'pyfrag' shell script is accessible via the terminal.:"
echo "$PYFRAGROOT/bin"