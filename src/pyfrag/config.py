"""
Configuration module for PyFrag package.

This module handles environment configuration, virtual environment setup,
and package-wide settings.
"""

import os
import sys
from pathlib import Path
from typing import Dict, Optional


class PyFragConfig:
    """Configuration manager for PyFrag."""

    def __init__(self):
        self._pyfrag_home = None
        self._virtual_env = None

    @property
    def pyfrag_home(self) -> Path:
        """Get the PyFrag home directory."""
        if self._pyfrag_home is None:
            # Try environment variable first
            if "PYFRAGHOME" in os.environ:
                self._pyfrag_home = Path(os.environ["PYFRAGHOME"])
            else:
                # Default to package location
                self._pyfrag_home = Path(__file__).parent.parent.parent.absolute()
        return self._pyfrag_home

    @pyfrag_home.setter
    def pyfrag_home(self, path: Path):
        """Set the PyFrag home directory."""
        self._pyfrag_home = Path(path)
        os.environ["PYFRAGHOME"] = str(self._pyfrag_home)

    @property
    def virtual_env(self) -> Optional[Path]:
        """Get the virtual environment path if available."""
        if self._virtual_env is None:
            # Check for virtual environment in the PyFrag directory
            venv_path = self.pyfrag_home / ".venv"
            if venv_path.exists():
                self._virtual_env = venv_path
            elif "VIRTUAL_ENV" in os.environ:
                self._virtual_env = Path(os.environ["VIRTUAL_ENV"])
        return self._virtual_env

    def setup_environment(self) -> Dict[str, str]:
        """Set up the environment variables for PyFrag."""
        env = os.environ.copy()
        env["PYFRAGHOME"] = str(self.pyfrag_home)
        env["HOSTPYFRAG"] = str(self.pyfrag_home)  # For compatibility with shell scripts

        # Add PyFrag scripts directory to PATH if not already there
        scripts_dir = str(self.get_scripts_path())
        current_path = env.get("PATH", "")
        if scripts_dir not in current_path:
            env["PATH"] = f"{scripts_dir}:{current_path}"

        # Add src directory to PYTHONPATH for imports
        src_dir = str(self.pyfrag_home / "src")
        current_pythonpath = env.get("PYTHONPATH", "")
        if src_dir not in current_pythonpath:
            if current_pythonpath:
                env["PYTHONPATH"] = f"{src_dir}:{current_pythonpath}"
            else:
                env["PYTHONPATH"] = src_dir

        return env

    def activate_virtual_env(self):
        """Activate the virtual environment if available."""
        if self.virtual_env:
            # Add virtual environment to Python path
            venv_site_packages = self.virtual_env / "lib" / f"python{sys.version_info.major}.{sys.version_info.minor}" / "site-packages"
            if venv_site_packages.exists() and str(venv_site_packages) not in sys.path:
                sys.path.insert(0, str(venv_site_packages))

            # Set virtual environment variables
            os.environ["VIRTUAL_ENV"] = str(self.virtual_env)
            os.environ["PATH"] = f"{self.virtual_env / 'bin'}:{os.environ.get('PATH', '')}"

    def get_adf_new_path(self) -> Path:
        """Get the path to the adf_new module."""
        return self.pyfrag_home / "src" / "pyfrag" / "host" / "standalone" / "adf_new"

    def get_host_path(self) -> Path:
        """Get the path to the host module."""
        return self.pyfrag_home / "src" / "pyfrag" / "host"

    def get_executables_path(self) -> Path:
        """Get the path to the executables (former host/bin)."""
        return self.get_host_path() / "bin"

    def get_scripts_path(self) -> Path:
        """Get the path to the main scripts."""
        return self.get_host_path() / "bin"

    def __str__(self) -> str:
        return f"PyFragConfig(pyfrag_home={self.pyfrag_home}, virtual_env={self.virtual_env})"


# Global configuration instance
config = PyFragConfig()


def get_config() -> PyFragConfig:
    """Get the global PyFrag configuration instance."""
    return config


def setup_pyfrag_environment():
    """Set up the PyFrag environment."""
    config.activate_virtual_env()
    return config.setup_environment()
