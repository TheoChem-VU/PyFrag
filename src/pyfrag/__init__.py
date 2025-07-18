"""
PyFrag: Automation of fragment calculations in quantum chemistry.

This package provides tools for automating fragment calculations using
various quantum chemistry engines including ADF, Gaussian, ORCA, and Turbomole.
"""

__version__ = "2025.1.0"
__author__ = "X. Sun, J. Poater, T. A. Hamlin, F. M. Bickelhaupt, S.J. Lekanne Deprez"

from . import config, executables
from .main import main

__all__ = ["main", "config", "executables", "__version__", "__author__"]
