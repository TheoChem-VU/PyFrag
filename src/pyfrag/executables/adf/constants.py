from typing import List

from scm.plams import Units

# =====================================================================
# Constants and settings
# =====================================================================

HA_TO_EV: float = Units.conversion_ratio("Ha", "eV")
HA_TO_KCAL: float = Units.conversion_ratio("Ha", "kcal/mol")
BOHR_TO_ANGSTROM: float = Units.conversion_ratio("Bohr", "Angstrom")
RAD_TO_DEG: float = Units.conversion_ratio("rad", "deg")

SYSTEM_NAMES: List[str] = ["complex", "frag1", "frag2"]
