#!/bin/bash

# ==================================================================
# This script is used to deploy the PyFrag application on a High performance computer (HPC) cluster.
# The following is required for deploying the PyFrag application:
#
#   1. Install the Python package manager "uv" if it is not already installed
#   2. Create a virtual environment (Python) with the pyfrag package installed through `uv sync`
#
# ==================================================================

echo "Trying to deploy PyFrag on the cluster which involves:"
echo "1. Installing the Python package manager 'uv'"
echo "2. Creating a virtual environment (Python) with the pyfrag package installed through 'uv sync'"
echo "3. Uninstalling uv"

# NOTE: this function assumes the Linux/MacOS operating system
# We do not include yet support for Windows.
install_uv() {
    # Try first with the curl command
    curl -LsSf https://astral.sh/uv/install.sh | sh

    # Otherwise, try with wget
    wget -qO- https://astral.sh/uv/install.sh | sh

    echo "uv installed successfully"
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
uv sync --no-cache --no-dev

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
rm -r "$(uv python dir)"
rm -r "$(uv tool dir)"
rm ~/.local/bin/uv ~/.local/bin/uvx

echo "Finally, please make sure to set the PYFRAGSOURCE environment variable in the module file to /src/pyfrag folder so that the `pyfrag` shell script is accessible via the terminal."