import logging
import re
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Sequence, Tuple, Union

import constants as const
import numpy as np
from input import InputKeys
from scm.plams import AMSJob, AMSResults, Molecule, Units, dihedral
from scm.plams.mol.atom import Atom

try:
    from .input import Angle, BondLength, Dihedral
except ImportError:
    from input import Angle, BondLength, Dihedral

logger = logging.getLogger("PyFragResults-ADF")

if TYPE_CHECKING:
    from input import OrbitalDescription


# ======================================================================
# Helper functions
# ======================================================================


def convert_output_data_into_seperate_keys(data: Dict[str, Any]):
    """
    This function takes a dictionary and converts it into a format where each value in a list
    is associated with a unique key.

    For example, it converts {'overlap': [1.55, 1.99]} into
    {'overlap_1': 1.55, 'overlap_2':1.99}.
    """
    output_table = {}

    for key, value in data.items():
        # If the value is a list with more than one item, create a new key for each item
        if isinstance(value, list) and len(value) > 1:
            for index, item in enumerate(value, start=1):
                new_key = f"{key}{index}"
                output_table[new_key] = item
        # If the value is a list with one item or not a list, keep the original key
        elif isinstance(value, list) and len(value) == 1 and key in ("bondlength", "angle"):
            new_key = f"{key}1"
            output_table[new_key] = value[0]

        else:
            output_table[key] = value[0] if isinstance(value, list) else value

    return output_table


def get_eda_terms(complex_result: AMSResults, unit: str = "kcal/mol") -> Dict[str, float]:
    """
    This function reads the (canonical) EDA terms from the complex_result and returns them as a dictionary.
    These include: Pauli, Electrostatic, Orbital Interaction, Dispersion, Bond Energy, and Total Energy.

    Args:
        complex_result (AMSResults): The AMS results object.
    Returns:
        Dict[str, float]: A dictionary containing the EDA terms.
    """
    name_to_variable_mapping = {
        "Pauli": "Pauli Total",
        "Elstat": "Electrostatic Interaction",
        "OI": "Orb.Int. Total",
        "Disp": "Dispersion Energy",
        "Int": "Bond Energy",
    }
    eda_terms = {}
    for name, variable in name_to_variable_mapping.items():
        try:
            value = float(complex_result.readrkf("Energy", variable, file="adf")) * Units.conversion_ratio("hartree", unit)  # type: ignore
            eda_terms[name] = value
        except KeyError:
            logger.warning(f"Key {variable} not found in the results when reading the EDA terms.")
            eda_terms[name] = np.nan
    return eda_terms


def ensure_list(obj) -> List[Any]:
    # single number in adf t21 is number fommat which is need to convert list
    return obj if isinstance(obj, list) else [obj]


def match_pattern(patterns: List[str], input_string: str, group_replacements: List[Tuple[str, str]]) -> Dict[str, Union[str, int]]:
    """
    Matches the input string against a list of patterns and extracts components based on group replacements.

    Args:
        patterns (List[str]): List of regex patterns to match.
        input_string (str): The string to match against the patterns.
        group_replacements (List[Tuple[str, str]]): List of tuples specifying regex replacements for groups.

    Returns:
        Dict[str, Union[str, int]]: A dictionary of extracted components.
    """
    for pattern in patterns:
        match_obj = re.match(pattern, input_string)
        if match_obj:
            result = {}
            for group_name, replacement in group_replacements:
                result[group_name] = re.sub(replacement, "", match_obj.group())
            return result
    return {group_name: "UNKNOWN" for group_name, _ in group_replacements}


def check_for_irrep_oi_keys(irreps_raw: Sequence[str]) -> List[Dict[str, str]]:
    """Checks whether OI irreps are present that can be included. Returns a list of dictionaries with irreps that can be included."""
    # Split degenerate irreps (e.g., E1:1, E1:2) into one entry per irrep (e.g., E1) as these are stored in the kf file
    irreps = set([irrep if ":" not in irrep else irrep.split(":")[0] for irrep in irreps_raw])
    irreps = sorted(irreps)
    logger.info(msg=f"Found irreps {', '.join(irreps)} in complex that will be included in OI")
    irreps = [{"irrep": irrep} for irrep in irreps]

    return irreps


def get_orbital_interaction_energy(complex_result: AMSResults, irreps_raw: Sequence[str], irrep: str, orb_int_term: float) -> float:
    irreps = list(set([irrep if ":" not in irrep else irrep.split(":")[0] for irrep in irreps_raw]))
    irrepOI: list[float] = [complex_result.readrkf("Energy", "Orb.Int. " + irrep, file="adf") * const.HA_TO_KCAL for irrep in irreps]  # type: ignore # No clear typing from plams

    total = sum(irrepOI)
    if np.isclose(total, 0):
        logger.debug("Sum of irrepOI is zero; cannot compute fitCoefficient, defaulting to 0.0.")
        return 0.0
    fitCoefficient = orb_int_term / total
    return fitCoefficient * irrepOI[irreps.index(irrep)]


def get_vdd_charges(complex_result: AMSResults, atomList: Sequence[int]) -> List[float]:
    vddList = []
    for atom in atomList:
        vddScf: float = ensure_list(complex_result.readrkf("Properties", "AtomCharge_SCF Voronoi", file="adf"))[int(atom) - 1]
        vddInit: float = ensure_list(complex_result.readrkf("Properties", "AtomCharge_initial Voronoi", file="adf"))[int(atom) - 1]
        vddList.append((vddScf - vddInit) * 1000)
    return vddList


def get_hirshfeld_charges(complex_result: AMSResults, fragment_index: int) -> List[float]:
    charges = ensure_list(complex_result.readrkf("Properties", "FragmentCharge Hirshfeld", file="adf"))
    return charges[fragment_index - 1] if isinstance(fragment_index, int) else np.nan


def get_fragment_orbital_irrep(irreps_raw: Sequence[str], irrep_orb_number: Sequence[int]) -> List[str]:
    # append complex irrep label to each orbital, if symmetry is A, convert self.irrepOrbNum which is float type into list
    fa_irrepone = [[irrep for _ in range(number)] for irrep, number in zip(irreps_raw, irrep_orb_number)]
    return [irrep for sublist in fa_irrepone for irrep in sublist]


def GetOrbNum(irrep_orb_number: Sequence[int], core_orb_number: Sequence[int]) -> List[int]:
    """
    Returns the global orbital indices for all orbitals (including frozen core orbitals) in the complex, across all irreps.

    This function is used to map the full list of orbital indices (1-based) for the entire complex, including both core and valence orbitals, as required for population analysis and other properties that refer to the full set of orbitals.

    Example output: [1,2,3,4,5,6,7,8,9,...] (all orbitals, core+valence, in order by irrep)

    Args:
        irrep_orb_number (Sequence[int]): Number of valence orbitals per irrep.
        core_orb_number (Sequence[int]): Number of core orbitals per irrep.

    Returns:
        List[int]: List of global orbital indices (1-based) for all orbitals in the complex.
    """
    orb_numbers = []
    orbSum = 0
    for nrShell, nrCore in zip(irrep_orb_number, core_orb_number):
        orbSum += nrShell + nrCore
        orb_numbers.extend(list(range(orbSum - nrShell + 1, orbSum + 1)))
    return orb_numbers


def GetFragOrbNum(irrep_orb_number: Sequence[int], core_orb_number: Sequence[int]) -> List[int]:
    """
    Returns the orbital indices for valence (non-core) orbitals only, skipping frozen core orbitals, for all irreps.

    This function is used to map the indices of fragment orbitals (SFOs) that are relevant for fragment-based properties, such as orbital overlaps and fragment populations. The indices start after the core orbitals for each irrep.

    Example output: [core+1, core+2, ..., core+valence] for each irrep, concatenated.

    Args:
        irrep_orb_number (Sequence[int]): Number of valence orbitals per irrep.
        core_orb_number (Sequence[int]): Number of core orbitals per irrep.

    Returns:
        List[int]: List of orbital indices (1-based) for valence (non-core) orbitals only, for all irreps.
    """
    orb_numbers = []
    for nrShell, nrCore in zip(irrep_orb_number, core_orb_number):
        orb_numbers.extend(range(nrCore + 1, nrShell + nrCore + 1))
    return orb_numbers


def selected_atom_indices_in_ordered_complex(atom_indices_per_fragment: Mapping[str, Sequence[int]], atom_indices: Sequence[int]) -> List[int]:
    """Change atom number in old presentation of a molecule into atom number in new presentation that formed by assembling fragments"""
    original_fragment_indices = [atomNum for key in sorted(list(atom_indices_per_fragment.keys())) for atomNum in list(atom_indices_per_fragment[key])]
    selected_indices = [i for i in atom_indices if i in original_fragment_indices]
    return selected_indices


def get_fragment_number(complex_result: AMSResults, frag: str) -> int:
    # change frag type like 'frag1' into number like "1" recorded in t21
    # append fragmenttype(like 1 or 2) to each orbital
    fragType = str(complex_result.readrkf("Geometry", "fragmenttype", file="adf")).split()
    return fragType.index(frag) + 1


def split_homo_lumo_index(orbSign: str) -> Dict[str, Union[str, int]]:
    """
    Converts HOMO/LUMO/HOMO-1/LUMO+1/INDEX into a dictionary {'holu': 'HOMO', 'num': -1}.
    """
    patterns = [r"HOMO(.*)", r"LUMO(.*)", r"INDEX"]
    group_replacements = [("holu", r"(.[0-9]+)"), ("num", r"([a-zA-Z]+)")]
    result = match_pattern(patterns, orbSign, group_replacements)
    result["num"] = int(result["num"]) if str(result["num"]).isdigit() else 0
    return result


# -----------------------------------------------------------------------------
# Shared results printing functions for obtaining the outputData dictionary.
# This dictionary contains all the results obtained from the inputKeys and
# will be printed in an output table
# -----------------------------------------------------------------------------


def get_bondlength_results(bond_lengths: Sequence[BondLength], fragment_indices: Mapping[str, Sequence[int]], complexMolecule: Molecule) -> Dict[str, float]:
    """Get bond length results for the specified bond lengths.

    Arguments:
        bond_lengths (Sequence[BondLength]): The bond lengths to process.
        fragment_indices (Mapping[str, Sequence[int]]): The mapping of fragment indices.
        complexMolecule (Molecule): The complex molecule to reference.

    Returns:
        Dict[str, float]: The label and its corresponding bond length result (in angstroms).

    Example:
        {
            "bondlength1_1C-2C": 0.1,
            "bondlength2_1C-3C": -0.2
        }
    """
    outputData = {}
    for i, bond_length in enumerate(bond_lengths, start=1):
        atom_indices = selected_atom_indices_in_ordered_complex(fragment_indices, bond_length["atom_indices"])
        atoms: List[Atom] = [complexMolecule[atom_index] for atom_index in atom_indices]
        label = f"bondlength{i}_{'-'.join([f'{str(atom_index)}{atom.symbol}' for atom_index, atom in zip(atom_indices, atoms)])}"
        outputData[label] = atoms[0].distance_to(atoms[1]) - bond_length["original_value"]
    return outputData


def get_angle_results(angles: Sequence[Angle], fragment_indices: Mapping[str, Sequence[int]], complexMolecule: Molecule) -> Dict[str, float]:
    """Get angle results for the specified angles.

    Arguments:
        angles (Sequence[Angle]): The angles to process.
        fragment_indices (Mapping[str, Sequence[int]]): The mapping of fragment indices.
        complexMolecule (Molecule): The complex molecule to reference.

    Returns:
        Dict[str, float]: The label and its corresponding angle result (in degrees).

    Important: the `atom.angle()` method takes the middle atom as the reference

    Example:
        {
            "angle1_1C-2C-3C": 5.0,
            "angle2_1C-2C-3C": -3.0
        }
    """
    outputData = {}
    for i, angle in enumerate(angles, start=1):
        atom_indices = selected_atom_indices_in_ordered_complex(fragment_indices, angle["atom_indices"])
        atoms: List[Atom] = [complexMolecule[atom_indices[0]], complexMolecule[atom_indices[1]], complexMolecule[atom_indices[2]]]
        label = f"angle{i}_{'-'.join([f'{str(atom_index)}{atom.symbol}' for atom_index, atom in zip(atom_indices, atoms)])}"
        outputData[label] = atoms[1].angle(atoms[0], atoms[2], result_unit="degree") - angle["original_value"]
    return outputData


def get_dihedral_results(dihedrals: Sequence[Dihedral], fragment_indices: Mapping[str, Sequence[int]], complexMolecule: Molecule) -> Dict[str, float]:
    """
    Get dihedral results for the specified dihedrals.

    Arguments:
        dihedrals (Sequence[Dihedral]): The dihedrals to process.
        fragment_indices (Mapping[str, Sequence[int]]): The mapping of fragment indices.
        complexMolecule (Molecule): The complex molecule to reference.

    Returns:
        Dict[str, float]: The label and its corresponding dihedral result (in degrees).

    Example:
        {
            "dihedral1_1C-2C-3C-4C": 5.0,
            "dihedral2_1C-2C-3C-4C": -3.0
        }
    """
    outputData = {}
    for i, dihedral_angle in enumerate(dihedrals, start=1):
        atom_indices = selected_atom_indices_in_ordered_complex(fragment_indices, dihedral_angle["atom_indices"])
        atoms: List[Atom] = [complexMolecule[atom_indices[0]], complexMolecule[atom_indices[1]], complexMolecule[atom_indices[2]], complexMolecule[atom_indices[3]]]
        label = f"dihedral{i}_{'-'.join([f'{str(atom_index)}{atom.symbol}' for atom_index, atom in zip(atom_indices, atoms)])}"
        outputData[label] = dihedral(atoms[0].coords, atoms[1].coords, atoms[2].coords, atoms[3].coords, unit="degree") - dihedral_angle["original_value"]
    return outputData


def get_fragment_orbital_energy(complex_result: AMSResults, orb_fragment: Sequence[int], frag_irrep: Sequence[str], frag_orb: Sequence[int], index: int) -> float:
    energy = np.atleast_1d(complex_result.readrkf(f"Ftyp {str(orb_fragment[index])}{frag_irrep[index]}", "eps", file="adf"))[frag_orb[index] - 1]
    return energy * const.HA_TO_EV


def get_vdd_output_results(fragment_index_mapping, vdd_indices: Sequence[int], complex_mol: Molecule, complex_result: AMSResults) -> Dict[str, float]:
    """Get VDD output results for each fragment.

    Returns:
        Dict[str, float]: VDD charges for each fragment.
    """
    outputData = {}
    atom_indices = selected_atom_indices_in_ordered_complex(fragment_index_mapping, vdd_indices)
    atoms = [complex_mol[atom_index] for atom_index in atom_indices]
    vdd_charges = get_vdd_charges(complex_result, atom_indices)
    for i, charge in enumerate(vdd_charges, start=1):
        outputData[f"VDD{i}_{atoms[i - 1].symbol}{atom_indices[i - 1]}"] = charge
    return outputData


# ======================================================================
# Functions for making print labels for the output table
# ======================================================================


def make_orbital_label(orbital: "OrbitalDescription", include_frag: bool = False) -> str:
    """
    Makes a label for the orbital based on its type, fragment, irrep, and index.
    Current formats are:
        (with HOMO/LUMO specification): "{type}"
        (with irrep/index specification): "{index} {irrep}"

    The `input_frag` argument is used to distinguish between `overlaps` (always printed in the order fragment 1 - fragment 2), and other types of orbital labels such as `orbitalenergy` and `population`.

    """
    label: List[str] = []

    if include_frag:
        label.append(orbital["frag"])

    # HOMO/LUMO specification
    if orbital["type"].upper() != "INDEX":
        label.append(orbital["type"].upper())
        return "-".join(label)

    # Irrep/index specification
    if orbital["index"] is not None and orbital["irrep"] is not None:
        irrep: str = orbital["irrep"]
        index: str = orbital["index"]
        spin: Union[str, None] = None

        # If the index contains an underscore, it means the spin is included
        # For example, "2_A"
        if "_" in index:
            index, spin = index.split("_")

        if spin is not None:
            irrep = f"{irrep}-{spin}"

        label.append(f"{index}{irrep}")

        return "-".join(label)

    raise ValueError(f"Invalid orbital description: {orbital}")


# ======================================================================
# Classes for reading ADF results and extracting data
# ======================================================================


class PyFragRestrictedResult:
    def __init__(self, complex_result: AMSResults, inputKeys: InputKeys):  # __init__(self, complexJob, inputKeys)
        # 1. needed for output requested by user 2. complexJob.check passes
        self.complex_result = complex_result
        # irrep label for symmetry of complex
        self.irreps_raw = str(complex_result.readrkf("Symmetry", "symlab", file="adf")).split()
        # orbital numbers according to the symmetry of each fragment and the orbitals belonging to the same symmetry in different fragments
        self.fragOrb: list[int] = ensure_list(complex_result.readrkf("SFOs", "ifo", file="adf"))
        # symmetry for each orbital of fragments
        self.fragIrrep: list[str] = ensure_list(str(complex_result.readrkf("SFOs", "subspecies", file="adf")).split())
        # the fragment label for each orbital
        self.orb_fragment: list[int] = ensure_list(complex_result.readrkf("SFOs", "fragment", file="adf"))
        # number of orbitals for each symmetry for complex
        self.irrep_orb_number = ensure_list(complex_result.readrkf("Symmetry", "norb", file="adf"))
        # number of core orbitals for each symmetry for complex
        self.core_orb_number = ensure_list(complex_result.readrkf("Symmetry", "ncbs", file="adf"))

        if "orbitalenergy" in inputKeys:
            try:
                self.orb_energy: list[float] = ensure_list(complex_result.readrkf("SFOs", "escale", file="adf"))
            except KeyError:
                logger.debug(msg="Reading non-relativistic orbital energies")
                self.orb_energy = ensure_list(complex_result.readrkf("SFOs", "energy", file="adf"))

        if "population" in inputKeys:
            # occupation of each orbitals which is either 0 or 2
            self.orb_occupation = ensure_list(complex_result.readrkf("SFOs", "occupation", file="adf"))

    def get_fragment_orbital_index(self, orbDescriptor: "OrbitalDescription"):
        # orbDescriptor = {'type' = "HOMO/LUMO/INDEX", 'frag'='#frag', 'irrep'='irrepname', 'index'=i}
        fragOrbnum: int = get_fragment_number(self.complex_result, orbDescriptor["frag"])  # get fragment number
        orb_index = 0
        if split_homo_lumo_index(orbDescriptor["type"])["holu"] == "HOMO":
            orb_index = sorted(range(len(self.orb_energy)), key=lambda x: self.orb_energy[x] if (self.orb_fragment[x] == fragOrbnum) and self.orb_occupation[x] != 0 else -1.0e100, reverse=True)[
                -int(split_homo_lumo_index(orbDescriptor["type"])["num"])
            ]
        elif split_homo_lumo_index(orbDescriptor["type"])["holu"] == "LUMO":
            orb_index = sorted(range(len(self.orb_energy)), key=lambda x: self.orb_energy[x] if (self.orb_fragment[x] == fragOrbnum) and self.orb_occupation[x] == 0 else +1.0e100)[int(split_homo_lumo_index(orbDescriptor["type"])["num"])]
        elif split_homo_lumo_index(orbDescriptor["type"])["holu"] == "INDEX" and orbDescriptor["index"] is not None:
            for i in range(len(self.orb_energy)):
                if self.orb_fragment[i] == fragOrbnum and self.fragIrrep[i] == orbDescriptor["irrep"] and self.fragOrb[i] == int(orbDescriptor["index"]):
                    orb_index: int = i
                    break
        return orb_index

    def get_sfo_overlap(self, index_1: int, index_2: int) -> float:
        # orbital numbers according to the symmetry of the complex
        fa_orb = GetFragOrbNum(self.irrep_orb_number, self.core_orb_number)
        fa_irrep = get_fragment_orbital_irrep(self.irreps_raw, self.irrep_orb_number)
        max_index = max(fa_orb[index_1], fa_orb[index_2])
        min_index = min(fa_orb[index_1], fa_orb[index_2])
        index = max_index * (max_index - 1) / 2 + min_index - 1
        if fa_irrep[index_1] == fa_irrep[index_2]:
            self.overlap_matrix: List[float] = ensure_list(self.complex_result.readrkf(fa_irrep[index_1], "S-CoreSFO", file="adf"))
            return abs(self.overlap_matrix[int(index)])
        else:
            return 0

    def read_sfo_population(self, index: int) -> float:
        orb_numbers: List[int] = GetOrbNum(self.irrep_orb_number, self.core_orb_number)
        # populations of all orbitals
        sfo_populs: List[float] = ensure_list(self.complex_result.readrkf("SFO popul", "sfo_grosspop", file="adf"))
        return sfo_populs[orb_numbers[index] - 1]

    def get_output_data(self, complexMolecule: Molecule, outputData: Dict[str, Any], inputKeys: InputKeys):
        # collect default energy parts for activation strain analysis
        outputData.update(get_eda_terms(self.complex_result))
        outputData["EnergyTotal"] = outputData["Int"] + outputData["StrainTotal"]

        # check for unspecified options such as irrep printing if not specified by user
        if not inputKeys["irrepOI"] and len(self.irreps_raw) != 1:
            inputKeys["irrepOI"] = check_for_irrep_oi_keys(self.irreps_raw)

        # collect user defined data
        if inputKeys["overlap"]:
            for i, od in enumerate(inputKeys["overlap"], start=1):
                od1 = self.get_fragment_orbital_index(od[0])
                od2 = self.get_fragment_orbital_index(od[1])
                label = f"overlap{i}_{make_orbital_label(od[0])}_{make_orbital_label(od[1])}"
                print(od1, od2, label)
                outputData[label] = self.get_sfo_overlap(od1, od2)

        if inputKeys["population"]:
            for i, od in enumerate(inputKeys["population"], start=1):
                label = f"population{i}_{make_orbital_label(od, include_frag=True)}"
                outputData[label] = self.read_sfo_population(self.get_fragment_orbital_index(od))

        if inputKeys["orbitalenergy"]:
            for i, od in enumerate(inputKeys["orbitalenergy"], start=1):
                label = f"orbitalenergy{i}_{make_orbital_label(od, include_frag=True)}"
                outputData[label] = get_fragment_orbital_energy(self.complex_result, self.orb_fragment, self.fragIrrep, self.fragOrb, self.get_fragment_orbital_index(od))

        if inputKeys["irrepOI"]:
            for od in inputKeys["irrepOI"]:
                outputData[f"irrepOI_{od['irrep']}"] = get_orbital_interaction_energy(self.complex_result, self.irreps_raw, od["irrep"], outputData["OI"])

        if inputKeys["VDD"]:
            outputData.update(get_vdd_output_results(inputKeys["fragment_indices"], inputKeys["VDD"], complexMolecule, self.complex_result))

        if inputKeys["hirshfeld"]:
            outputData["hirshfeld"] = [get_hirshfeld_charges(self.complex_result, get_fragment_number(self.complex_result, fragment_specifier)) for fragment_specifier in inputKeys["hirshfeld"]]

        if inputKeys["bondlength"]:
            outputData.update(get_bondlength_results(inputKeys["bondlength"], inputKeys["fragment_indices"], complexMolecule))

        if inputKeys["angle"]:
            outputData.update(get_angle_results(inputKeys["angle"], inputKeys["fragment_indices"], complexMolecule))

        if inputKeys["dihedral"]:
            outputData.update(get_dihedral_results(inputKeys["dihedral"], inputKeys["fragment_indices"], complexMolecule))

        return convert_output_data_into_seperate_keys(outputData)


class PyFragUnrestrictedResult:
    def __init__(self, complex_result: AMSResults, inputKeys: InputKeys):  # __init__(self, complexJob, inputKeys)
        self.complex_result = complex_result
        # irrep label for symmetry of complex
        self.irreps_raw = str(complex_result.readrkf("Symmetry", "symlab", file="adf")).split()
        # number of orbitals for each symmetry for complex
        self.irrep_orb_number = ensure_list(complex_result.readrkf("Symmetry", "norb", file="adf"))
        # number of core orbitals for each symmetry for complex
        self.core_orb_number = ensure_list(complex_result.readrkf("Symmetry", "ncbs", file="adf"))
        # orbital numbers according to the symmetry of each fragment and the orbitals belonging to the same symmetry in different fragments
        self.fragOrb: list[int] = ensure_list(complex_result.readrkf("SFOs", "ifo", file="adf"))
        # symmetry for each orbital of fragments
        self.fragIrrep: list[str] = ensure_list(str(complex_result.readrkf("SFOs", "subspecies", file="adf")).split())
        # the fragment label for each orbital
        self.orb_fragment: list[int] = ensure_list(complex_result.readrkf("SFOs", "fragment", file="adf"))

        if "orbitalenergy" in inputKeys:
            # energy for each orbital of spin A and B (escale is only if relativistic corrections are used)
            try:
                logging.log(level=logging.DEBUG, msg="Reading relativistic orbital energies")
                self.orb_energy = ensure_list(complex_result.readrkf("SFOs", "escale", file="adf"))
                self.orb_energy_B = ensure_list(complex_result.readrkf("SFOs", "escale_B", file="adf"))
            except KeyError:
                logging.log(level=logging.DEBUG, msg="Reading non-relativistic orbital energies")
                self.orb_energy = ensure_list(complex_result.readrkf("SFOs", "energy", file="adf"))
                self.orb_energy_B = ensure_list(complex_result.readrkf("SFOs", "energy_B", file="adf"))

        if "population" in inputKeys:
            # occupation of each orbitals of A which is either 0 or 2
            self.orb_occupation = ensure_list(complex_result.readrkf("SFOs", "occupation", file="adf"))
            # occupation of each orbitals of b which is either 0 or 2
            self.orb_occupation_B = ensure_list(complex_result.readrkf("SFOs", "occupation_B", file="adf"))

    def split_orbital_into_spin_and_index(self, orbSign: str) -> Dict[str, str]:
        # convert 1_A or 1_B  into dict {'spin': 'A', 'num': -1}
        for matchString in [r"(.*)_A", r"(.*)_B"]:
            matchObj = re.match(matchString, orbSign)
            if matchObj:
                spin = re.sub(r"([0-9]+_)", "", matchObj.group())
                num = re.sub(r"(_[a-zA-Z]+)", "", matchObj.group())
                return {"spin": spin, "num": num}
        # Default return to ensure all code paths return a value
        return {"spin": "UNKNOWN", "num": "0"}

    def get_fragment_orbital_index(self, orbDescriptor: "OrbitalDescription") -> Tuple[int, str]:
        # orbDescriptor = {'type' = "HOMO/LUMO/INDEX", 'frag'='#frag', 'irrep'='irrepname', 'index'=i}
        fragOrbnum = get_fragment_number(self.complex_result, orbDescriptor["frag"])  # get fragment number

        orb_index = 0
        spin = ""
        orb_index_AB = 0
        orb_energy: List[float] = self.orb_energy + self.orb_energy_B
        orb_fragment: List[int] = self.orb_fragment + self.orb_fragment
        orb_occupation: List[float] = self.orb_occupation + self.orb_occupation_B

        # for spin A
        if split_homo_lumo_index(orbDescriptor["type"])["holu"] == "HOMO":
            orb_index_AB = sorted(range(len(self.orb_energy) * 2), key=lambda x: orb_energy[x] if (orb_fragment[x] == fragOrbnum) and orb_occupation[x] != 0 else -1.0e100, reverse=True)[-int(split_homo_lumo_index(orbDescriptor["type"])["num"])]

            if orb_index_AB <= len(self.orb_energy) - 1:
                orb_index = orb_index_AB
                spin = "A"
            else:
                orb_index = orb_index_AB - len(self.orb_energy)
                spin = "B"

        elif split_homo_lumo_index(orbDescriptor["type"])["holu"] == "LUMO":
            orb_index_AB = sorted(range(len(self.orb_energy) * 2), key=lambda x: orb_energy[x] if (orb_fragment[x] == fragOrbnum) and orb_occupation[x] == 0 else +1.0e100)[int(split_homo_lumo_index(orbDescriptor["type"])["num"])]

            if orb_index_AB <= len(self.orb_energy) - 1:
                orb_index = orb_index_AB
                spin = "A"
            else:
                orb_index = orb_index_AB - len(self.orb_energy)
                spin = "B"

        elif split_homo_lumo_index(orbDescriptor["type"])["holu"] == "INDEX":
            spinOrbnum = int(self.split_orbital_into_spin_and_index(str(orbDescriptor["index"]))["num"])
            spinOrbspin = self.split_orbital_into_spin_and_index(str(orbDescriptor["index"]))["spin"]
            if spinOrbspin == "A":
                for i in range(len(self.orb_energy)):
                    if self.orb_fragment[i] == fragOrbnum and self.fragIrrep[i] == orbDescriptor["irrep"] and self.fragOrb[i] == spinOrbnum:
                        orb_index = i
                        spin = "A"
                        break
            if spinOrbspin == "B":
                for i in range(len(self.orb_energy_B)):
                    if self.orb_fragment[i] == fragOrbnum and self.fragIrrep[i] == orbDescriptor["irrep"] and self.fragOrb[i] == spinOrbnum:
                        orb_index = i
                        spin = "B"
                        break

        return (orb_index, spin)

    def get_sfo_overlap(self, index_1: Tuple[int, str], index_2: Tuple[int, str]) -> float:
        # orbital numbers according to the symmetry of the complex
        fa_orb = GetFragOrbNum(self.irrep_orb_number, self.core_orb_number)
        fa_irrep = get_fragment_orbital_irrep(self.irreps_raw, self.irrep_orb_number)
        max_index = max(fa_orb[index_1[0]], fa_orb[index_2[0]])
        min_index = min(fa_orb[index_1[0]], fa_orb[index_2[0]])
        index = max_index * (max_index - 1) / 2 + min_index - 1
        if fa_irrep[index_1[0]] == fa_irrep[index_2[0]]:
            if index_1[1] == "A" and index_2[1] == "A":
                self.overlap_matrix: List[float] = ensure_list(self.complex_result.readrkf(fa_irrep[index_1[0]], "S-CoreSFO", file="adf"))
                return abs(self.overlap_matrix[int(index)])
            elif index_1[1] == "B" and index_2[1] == "B":
                self.overlap_matrix: List[float] = ensure_list(self.complex_result.readrkf(fa_irrep[index_1[0]], "S-CoreSFO_B", file="adf"))
                return abs(self.overlap_matrix[int(index)])
            else:
                return 0
        else:
            return 0

    def read_sfo_population(self, index: Tuple[int, str]) -> float:
        orb_numbers = GetOrbNum(self.irrep_orb_number, self.core_orb_number)
        # populations of all orbitals
        sfo_populs: List[float] = ensure_list(self.complex_result.readrkf("SFO popul", "sfo_grosspop", file="adf"))
        if index[1] == "A":
            return sfo_populs[orb_numbers[index[0]] - 1]
        else:
            return sfo_populs[orb_numbers[index[0]] + len(self.orb_energy) - 1]

    def get_output_data(self, complexMolecule: Molecule, outputData: Dict[str, Any], inputKeys: InputKeys):
        # collect default energy parts for activation strain analysis
        outputData.update(get_eda_terms(self.complex_result))
        outputData["EnergyTotal"] = outputData["Int"] + outputData["StrainTotal"]

        # check for unspecified options such as irrep printing if not specified by user
        if not inputKeys["irrepOI"] and len(self.irreps_raw) != 1:
            inputKeys["irrepOI"] = check_for_irrep_oi_keys(self.irreps_raw)

        # collect user defined data
        if inputKeys["overlap"]:
            for i, od in enumerate(inputKeys["overlap"], start=1):
                od1 = self.get_fragment_orbital_index(od[0])
                od2 = self.get_fragment_orbital_index(od[1])
                label = f"overlap{i}_{make_orbital_label(od[0])}_{make_orbital_label(od[1])}"
                outputData[label] = self.get_sfo_overlap(od1, od2)

        if inputKeys["population"]:
            for i, od in enumerate(inputKeys["population"], start=1):
                label = f"population{i}_{make_orbital_label(od, include_frag=True)}"
                outputData[label] = self.read_sfo_population(self.get_fragment_orbital_index(od))

        if inputKeys["orbitalenergy"]:
            for i, od in enumerate(inputKeys["orbitalenergy"], start=1):
                label = f"orbitalenergy{i}_{make_orbital_label(od, include_frag=True)}"
                outputData[label] = get_fragment_orbital_energy(self.complex_result, self.orb_fragment, self.fragIrrep, self.fragOrb, self.get_fragment_orbital_index(od)[0])

        if inputKeys["irrepOI"]:
            for od in inputKeys["irrepOI"]:
                outputData[f"irrepOI_{od['irrep']}"] = get_orbital_interaction_energy(self.complex_result, self.irreps_raw, od["irrep"], outputData["OI"])

        if inputKeys["VDD"]:
            outputData.update(get_vdd_output_results(inputKeys["fragment_indices"], inputKeys["VDD"], complexMolecule, self.complex_result))

        if inputKeys["hirshfeld"]:
            outputData["hirshfeld"] = [get_hirshfeld_charges(self.complex_result, get_fragment_number(self.complex_result, fragment_specifier)) for fragment_specifier in inputKeys["hirshfeld"]]

        if inputKeys["bondlength"]:
            outputData.update(get_bondlength_results(inputKeys["bondlength"], inputKeys["fragment_indices"], complexMolecule))

        if inputKeys["angle"]:
            outputData.update(get_angle_results(inputKeys["angle"], inputKeys["fragment_indices"], complexMolecule))

        if inputKeys["dihedral"]:
            outputData.update(get_dihedral_results(inputKeys["dihedral"], inputKeys["fragment_indices"], complexMolecule))

        return convert_output_data_into_seperate_keys(outputData)


# ======================================================================
# Interface functions
# ======================================================================


def get_pyfrag_results(complexJob: AMSJob, inputKeys: InputKeys) -> Union[PyFragRestrictedResult, "PyFragUnrestrictedResult"]:
    """
    Returns a results instance of the PyFragResult class.
    The function checks if the job is a restricted or unrestricted calculation and returns the appropriate class instance.

    Args:
        complexJob (AMSJob): The AMS job object.
        inputKeys (Dict[str, Any]): A dictionary containing input keys.
    Returns:
        Union[PyFragRestrictedResult, PyFragUnrestrictedResult]: An instance of the appropriate results class.

    """
    adf_settings = complexJob.settings.input.adf

    unrestricted = False
    if any(key in adf_settings for key in ["spinpolarization", "unrestricted"]):
        logger.info("Unrestricted calculation detected.")
        unrestricted = True

    return PyFragUnrestrictedResult(complexJob.results, inputKeys) if unrestricted else PyFragRestrictedResult(complexJob.results, inputKeys)
