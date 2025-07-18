import logging
import re
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Sequence, Tuple, Union

import constants as const
import numpy as np
from input import InputKeys
from scm.plams import AMSJob, AMSResults, Molecule, Units
from scm.plams.mol.atom import Atom

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


def get_eda_terms(complexResult: AMSResults, unit: str = "kcal/mol") -> Dict[str, float]:
    """
    This function reads the (canonical) EDA terms from the complexResult and returns them as a dictionary.
    These include: Pauli, Electrostatic, Orbital Interaction, Dispersion, Bond Energy, and Total Energy.

    Args:
        complexResult (AMSResults): The AMS results object.
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
            value = float(complexResult.readrkf("Energy", variable, file="adf")) * Units.conversion_ratio("hartree", unit)  # type: ignore
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


def check_for_irrep_oi_keys(irrepType: Sequence[str]) -> List[Dict[str, str]]:
    """Checks whether OI irreps are present that can be included. Returns a list of dictionaries with irreps that can be included."""
    # Split degenerate irreps (e.g., E1:1, E1:2) into one entry per irrep (e.g., E1) as these are stored in the kf file
    irreps = set([irrep if ":" not in irrep else irrep.split(":")[0] for irrep in irrepType])
    irreps = sorted(irreps)
    logger.info(msg=f"Found irreps {', '.join(irreps)} in complex that will be included in OI")
    irreps = [{"irrep": irrep} for irrep in irreps]

    return irreps


def get_orbital_interaction_energy(complex_result: AMSResults, irrepType: Sequence[str], irrep: str, orb_int_term: float) -> float:
    irreps = list(set([irrep if ":" not in irrep else irrep.split(":")[0] for irrep in irrepType]))
    irrepOI: list[float] = [complex_result.readrkf("Energy", "Orb.Int. " + irrep, file="adf") * const.HA_TO_KCAL for irrep in irreps]  # type: ignore # No clear typing from plams
    fitCoefficient = orb_int_term / sum(irrepOI)  # MetaGGAs have a scaling factor that needs to be applied
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


def get_fragment_orbital_irrep(irrepType: Sequence[str], irrep_orb_number: Sequence[int]) -> List[str]:
    # append complex irrep label to each orbital, if symmetry is A, convert self.irrepOrbNum which is float type into list
    faIrrepone = [[irrep for _ in range(number)] for irrep, number in zip(irrepType, irrep_orb_number)]
    return [irrep for sublist in faIrrepone for irrep in sublist]


def GetOrbNum(irrep_orb_number: Sequence[int], core_orb_number: Sequence[int]) -> List[int]:
    # GetOrbNumbers including frozen core orbitals, this is necessary to read population, list like [3,4,5,12,13,14,15,23,24,26]
    # core orbital number corresponding to each irrep of complex symmetry
    orbNumbers = []
    orbSum = 0
    for nrShell, nrCore in zip(irrep_orb_number, core_orb_number):
        orbSum += nrShell + nrCore
        orbNumbers.extend(list(range(orbSum - nrShell + 1, orbSum + 1)))
    return orbNumbers


def GetFragOrbNum(irrep_orb_number: Sequence[int], core_orb_number: Sequence[int]) -> List[int]:
    # GetOrbNumbers including frozen core orbitals, this is necessary to read population, list like [3,4,5,12,13,14,15,23,24,26]
    # core orbital number corresponding to each irrep of complex symmetry
    orbNumbers = []
    for nrShell, nrCore in zip(irrep_orb_number, core_orb_number):
        orbNumbers.extend(range(nrCore + 1, nrShell + nrCore + 1))
    return orbNumbers


def get_atom_indices(atom_indices_per_fragment: Mapping[str, Sequence[int]], atom_indices: Sequence[int]) -> List[int]:
    # change atom number in old presentation of a molecule into atom number in new presentation that formed by assembling fragments
    original_fragment_indices = [atomNum for key in sorted(list(atom_indices_per_fragment.keys())) for atomNum in list(atom_indices_per_fragment[key])]
    reordered_indices = [original_fragment_indices.index(i) + 1 for i in atom_indices]
    return reordered_indices


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


# ======================================================================
# Functions for making print labels for the output table
# ======================================================================


def make_orbital_label(orbital: "OrbitalDescription", include_frag: bool = False) -> str:
    """
    Makes a label for the orbital based on its type, fragment, irrep, and index.
    Current formats are:
        (with HOMO/LUMO specification): "{type}"
        (with irrep/index specification): "{index} {irrep}"

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
        irrep = orbital["irrep"]

        # If the index contains an underscore, it means the spin is included
        if "_" in orbital["index"]:
            index, spin = orbital["index"].split("_")

            if spin is not None:
                irrep = f"{irrep}-{spin}"
            label.append(f"{index}{irrep}")

            return "-".join(label)

    return f"{orbital['frag']}"


# ======================================================================
# Classes for reading ADF results and extracting data
# ======================================================================


class PyFragRestrictedResult:
    def __init__(self, complexResult: AMSResults, inputKeys: InputKeys):  # __init__(self, complexJob, inputKeys)
        # 1. needed for output requested by user 2. complexJob.check passes
        self.complexResult = complexResult
        # irrep label for symmetry of complex
        self.irrepType = str(complexResult.readrkf("Symmetry", "symlab", file="adf")).split()

        for key in list(inputKeys.keys()):
            if key == "overlap" or key == "population" or key == "orbitalenergy" or key == "irrepOI":
                # orbital numbers according to the symmetry of each fragment and the orbitals belonging to the same symmetry in different fragments
                self.fragOrb: list[int] = ensure_list(complexResult.readrkf("SFOs", "ifo", file="adf"))
                # symmetry for each orbital of fragments
                self.fragIrrep: list[str] = ensure_list(str(complexResult.readrkf("SFOs", "subspecies", file="adf")).split())
                # the fragment label for each orbital
                self.orbFragment: list[int] = ensure_list(complexResult.readrkf("SFOs", "fragment", file="adf"))

                try:
                    self.orbEnergy: list[float] = ensure_list(complexResult.readrkf("SFOs", "escale", file="adf"))
                except KeyError:
                    logger.debug(msg="Reading non-relativistic orbital energies")
                    self.orbEnergy = ensure_list(complexResult.readrkf("SFOs", "energy", file="adf"))

                # occupation of each orbitals which is either 0 or 2
                self.orbOccupation = ensure_list(complexResult.readrkf("SFOs", "occupation", file="adf"))
                # number of orbitals for each symmetry for complex
                self.irrep_orb_number = ensure_list(complexResult.readrkf("Symmetry", "norb", file="adf"))
                # number of core orbitals for each symmetry for complex
                self.core_orb_number = ensure_list(complexResult.readrkf("Symmetry", "ncbs", file="adf"))

    def get_fragment_orbital_index(self, orbDescriptor: "OrbitalDescription"):
        # orbDescriptor = {'type' = "HOMO/LUMO/INDEX", 'frag'='#frag', 'irrep'='irrepname', 'index'=i}
        fragOrbnum: int = get_fragment_number(self.complexResult, orbDescriptor["frag"])  # get fragment number
        orbIndex = 0
        if split_homo_lumo_index(orbDescriptor["type"])["holu"] == "HOMO":
            orbIndex = sorted(range(len(self.orbEnergy)), key=lambda x: self.orbEnergy[x] if (self.orbFragment[x] == fragOrbnum) and self.orbOccupation[x] != 0 else -1.0e100, reverse=True)[-int(split_homo_lumo_index(orbDescriptor["type"])["num"])]
        elif split_homo_lumo_index(orbDescriptor["type"])["holu"] == "LUMO":
            orbIndex = sorted(range(len(self.orbEnergy)), key=lambda x: self.orbEnergy[x] if (self.orbFragment[x] == fragOrbnum) and self.orbOccupation[x] == 0 else +1.0e100)[int(split_homo_lumo_index(orbDescriptor["type"])["num"])]
        elif split_homo_lumo_index(orbDescriptor["type"])["holu"] == "INDEX" and orbDescriptor["index"] is not None:
            for i in range(len(self.orbEnergy)):
                if self.orbFragment[i] == fragOrbnum and self.fragIrrep[i] == orbDescriptor["irrep"] and self.fragOrb[i] == int(orbDescriptor["index"]):
                    orbIndex: int = i
                    break
        return orbIndex

    def get_sfo_overlap(self, index_1, index_2) -> float:
        # orbital numbers according to the symmetry of the complex
        faOrb = GetFragOrbNum(self.irrep_orb_number, self.core_orb_number)
        faIrrep = get_fragment_orbital_irrep(self.irrepType, self.irrep_orb_number)
        maxIndex = max(faOrb[index_1], faOrb[index_2])
        minIndex = min(faOrb[index_1], faOrb[index_2])
        index = maxIndex * (maxIndex - 1) / 2 + minIndex - 1
        if faIrrep[index_1] == faIrrep[index_2]:
            self.overlap_matrix: List[float] = ensure_list(self.complexResult.readrkf(faIrrep[index_1], "S-CoreSFO", file="adf"))
            return abs(self.overlap_matrix[int(index)])
        else:
            return 0

    def get_fragment_orbital_energy(self, index: int) -> float:
        return ensure_list(self.complexResult.readrkf("Ftyp " + str(self.orbFragment[index]) + self.fragIrrep[index], "eps", file="adf"))[self.fragOrb[index] - 1]

    def read_sfo_population(self, index) -> float:
        orbNumbers: List[int] = GetOrbNum(self.irrep_orb_number, self.core_orb_number)
        # populations of all orbitals
        sfoPopul: List[float] = ensure_list(self.complexResult.readrkf("SFO popul", "sfo_grosspop", file="adf"))
        return sfoPopul[orbNumbers[index] - 1]

    def get_output_data(self, complexMolecule: Molecule, outputData: Dict[str, Any], inputKeys: InputKeys):
        # collect default energy parts for activation strain analysis
        outputData.update(get_eda_terms(self.complexResult))
        outputData["EnergyTotal"] = outputData["Int"] + outputData["StrainTotal"]

        # check for unspecified options such as irrep printing if not specified by user
        if not inputKeys["irrepOI"] and len(self.irrepType) != 1:
            inputKeys["irrepOI"] = check_for_irrep_oi_keys(self.irrepType)

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
                outputData[label] = self.get_fragment_orbital_energy(self.get_fragment_orbital_index(od))

        if inputKeys["irrepOI"]:
            for od in inputKeys["irrepOI"]:
                outputData[f"irrepOI_{od['irrep']}"] = get_orbital_interaction_energy(self.complexResult, self.irrepType, od["irrep"], outputData["OI"])

        if inputKeys["VDD"]:
            outputData["VDD"] = get_vdd_charges(self.complexResult, get_atom_indices(inputKeys["fragment_indices"], inputKeys["VDD"]))

        if inputKeys["hirshfeld"]:
            outputData["hirshfeld"] = [get_hirshfeld_charges(self.complexResult, get_fragment_number(self.complexResult, fragment_specifier)) for fragment_specifier in inputKeys["hirshfeld"]]

        if inputKeys["bondlength"]:
            for i, od in enumerate(inputKeys["bondlength"], start=1):
                atom_indices = get_atom_indices(inputKeys["fragment_indices"], od["bond_definition"])
                atoms: List[Atom] = [complexMolecule[atom_indices[0]], complexMolecule[atom_indices[1]]]
                label = f"bondlength{i}_{'-'.join([f'{str(atom_index)}{atom.symbol}' for atom_index, atom in zip(atom_indices, atoms)])}"
                outputData[label] = atoms[0].distance_to(atoms[1]) - od["original_value"]

        if inputKeys["angle"]:
            for i, od in enumerate(inputKeys["angle"], start=1):
                atom_indices = get_atom_indices(inputKeys["fragment_indices"], od["angle_definition"])
                atoms: List[Atom] = [complexMolecule[atom_indices[0]], complexMolecule[atom_indices[1]], complexMolecule[atom_indices[2]]]
                label = f"angle{i}_{'-'.join([f'{str(atom_index)}{atom.symbol}' for atom_index, atom in zip(atom_indices, atoms)])}"
                outputData["angle"] = atoms[0].angle(atoms[1], atoms[2]) - od["original_value"]

        return convert_output_data_into_seperate_keys(outputData)


class PyFragUnrestrictedResult:
    def __init__(self, complexResult: AMSResults, inputKeys: InputKeys):  # __init__(self, complexJob, inputKeys)
        # 1. needed for output requested by user 2. complexJob.check passes
        self.complexResult = complexResult
        # irrep label for symmetry of complex
        self.irrepType = str(complexResult.readrkf("Symmetry", "symlab", file="adf")).split()

        for key in list(inputKeys.keys()):
            if key == "overlap" or key == "population" or key == "orbitalenergy" or key == "irrepOI":
                # orbital numbers according to the symmetry of each fragment and the orbitals belonging to the same symmetry in different fragments
                self.fragOrb: list[int] = ensure_list(complexResult.readrkf("SFOs", "ifo", file="adf"))
                # symmetry for each orbital of fragments
                self.fragIrrep: list[str] = ensure_list(str(complexResult.readrkf("SFOs", "subspecies", file="adf")).split())
                # the fragment label for each orbital
                self.orbFragment: list[int] = ensure_list(complexResult.readrkf("SFOs", "fragment", file="adf"))

                # energy for each orbital of spin A and B (escale is only if relativistic corrections are used)
                try:
                    logging.log(level=logging.DEBUG, msg="Reading relativistic orbital energies")
                    self.orbEnergy = ensure_list(complexResult.readrkf("SFOs", "escale", file="adf"))
                    self.orbEnergy_B = ensure_list(complexResult.readrkf("SFOs", "escale_B", file="adf"))
                except KeyError:
                    logging.log(level=logging.DEBUG, msg="Reading non-relativistic orbital energies")
                    self.orbEnergy = ensure_list(complexResult.readrkf("SFOs", "energy", file="adf"))
                    self.orbEnergy_B = ensure_list(complexResult.readrkf("SFOs", "energy_B", file="adf"))

                # occupation of each orbitals of A which is either 0 or 2
                self.orbOccupation = ensure_list(complexResult.readrkf("SFOs", "occupation", file="adf"))
                # occupation of each orbitals of b which is either 0 or 2
                self.orbOccupation_B = ensure_list(complexResult.readrkf("SFOs", "occupation_B", file="adf"))
                # number of orbitals for each symmetry for complex
                self.irrep_orb_number = ensure_list(complexResult.readrkf("Symmetry", "norb", file="adf"))
                # number of core orbitals for each symmetry for complex
                self.core_orb_number = ensure_list(complexResult.readrkf("Symmetry", "ncbs", file="adf"))

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
        fragOrbnum = get_fragment_number(self.complexResult, orbDescriptor["frag"])  # get fragment number

        orbIndex = 0
        spin = ""
        orbIndex_AB = 0
        orbEnergy = self.orbEnergy + self.orbEnergy_B
        orbFragment = self.orbFragment + self.orbFragment
        orbOccupation = self.orbOccupation + self.orbOccupation_B

        # for spin A
        if split_homo_lumo_index(orbDescriptor["type"])["holu"] == "HOMO":
            orbIndex_AB = sorted(range(len(self.orbEnergy) * 2), key=lambda x: orbEnergy[x] if (orbFragment[x] == fragOrbnum) and orbOccupation[x] != 0 else -1.0e100, reverse=True)[-int(split_homo_lumo_index(orbDescriptor["type"])["num"])]

            if orbIndex_AB <= len(self.orbEnergy) - 1:
                orbIndex = orbIndex_AB
                spin = "A"
            else:
                orbIndex = orbIndex_AB - len(self.orbEnergy)
                spin = "B"

        elif split_homo_lumo_index(orbDescriptor["type"])["holu"] == "LUMO":
            orbIndex_AB = sorted(range(len(self.orbEnergy) * 2), key=lambda x: orbEnergy[x] if (orbFragment[x] == fragOrbnum) and orbOccupation[x] == 0 else +1.0e100)[int(split_homo_lumo_index(orbDescriptor["type"])["num"])]

            if orbIndex_AB <= len(self.orbEnergy) - 1:
                orbIndex = orbIndex_AB
                spin = "A"
            else:
                orbIndex = orbIndex_AB - len(self.orbEnergy)
                spin = "B"

        elif split_homo_lumo_index(orbDescriptor["type"])["holu"] == "INDEX":
            spinOrbnum = int(self.split_orbital_into_spin_and_index(str(orbDescriptor["index"]))["num"])
            spinOrbspin = self.split_orbital_into_spin_and_index(str(orbDescriptor["index"]))["spin"]
            if spinOrbspin == "A":
                for i in range(len(self.orbEnergy)):
                    if self.orbFragment[i] == fragOrbnum and self.fragIrrep[i] == orbDescriptor["irrep"] and self.fragOrb[i] == spinOrbnum:
                        orbIndex = i
                        spin = "A"
                        break
            if spinOrbspin == "B":
                for i in range(len(self.orbEnergy_B)):
                    if self.orbFragment[i] == fragOrbnum and self.fragIrrep[i] == orbDescriptor["irrep"] and self.fragOrb[i] == spinOrbnum:
                        orbIndex = i
                        spin = "B"
                        break

        return (orbIndex, spin)

    def get_sfo_overlap(self, index_1: Tuple[int, str], index_2: Tuple[int, str]) -> float:
        # orbital numbers according to the symmetry of the complex
        faOrb = GetFragOrbNum(self.irrep_orb_number, self.core_orb_number)
        faIrrep = get_fragment_orbital_irrep(self.irrepType, self.irrep_orb_number)
        maxIndex = max(faOrb[index_1[0]], faOrb[index_2[0]])
        minIndex = min(faOrb[index_1[0]], faOrb[index_2[0]])
        index = maxIndex * (maxIndex - 1) / 2 + minIndex - 1
        if faIrrep[index_1[0]] == faIrrep[index_2[0]]:
            if index_1[1] == "A" and index_2[1] == "A":
                self.overlap_matrix: List[float] = ensure_list(self.complexResult.readrkf(faIrrep[index_1[0]], "S-CoreSFO", file="adf"))
                return abs(self.overlap_matrix[int(index)])
            elif index_1[1] == "B" and index_2[1] == "B":
                self.overlap_matrix: List[float] = ensure_list(self.complexResult.readrkf(faIrrep[index_1[0]], "S-CoreSFO_B", file="adf"))
                return abs(self.overlap_matrix[int(index)])
            else:
                return 0
        else:
            return 0

    def get_fragment_orbital_energy(self, index: int) -> float:
        energy = ensure_list(self.complexResult.readrkf("Ftyp " + str(self.orbFragment[index]) + self.fragIrrep[index], "eps", file="adf"))[self.fragOrb[index] - 1]
        return energy * const.HA_TO_EV

    def read_sfo_population(self, index: Tuple[int, str]) -> float:
        orbNumbers = GetOrbNum(self.irrep_orb_number, self.core_orb_number)
        # populations of all orbitals
        sfoPopul: List[float] = ensure_list(self.complexResult.readrkf("SFO popul", "sfo_grosspop", file="adf"))
        if index[1] == "A":
            return sfoPopul[orbNumbers[index[0]] - 1]
        else:
            return sfoPopul[orbNumbers[index[0]] + len(self.orbEnergy) - 1]

    def get_output_data(self, complexMolecule: Molecule, outputData: Dict[str, Any], inputKeys: InputKeys):
        # collect default energy parts for activation strain analysis
        outputData.update(get_eda_terms(self.complexResult))
        outputData["EnergyTotal"] = outputData["Int"] + outputData["StrainTotal"]

        # check for unspecified options such as irrep printing if not specified by user
        if not inputKeys["irrepOI"] and len(self.irrepType) != 1:
            inputKeys["irrepOI"] = check_for_irrep_oi_keys(self.irrepType)

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
                outputData[label] = self.get_fragment_orbital_energy(self.get_fragment_orbital_index(od)[0])

        if inputKeys["irrepOI"]:
            for od in inputKeys["irrepOI"]:
                outputData[f"irrepOI_{od['irrep']}"] = get_orbital_interaction_energy(self.complexResult, self.irrepType, od["irrep"], outputData["OI"])

        if inputKeys["VDD"]:
            outputData["VDD"] = get_vdd_charges(self.complexResult, get_atom_indices(inputKeys["fragment_indices"], inputKeys["VDD"]))

        if inputKeys["hirshfeld"]:
            outputData["hirshfeld"] = [get_hirshfeld_charges(self.complexResult, get_fragment_number(self.complexResult, fragment_specifier)) for fragment_specifier in inputKeys["hirshfeld"]]

        if inputKeys["bondlength"]:
            for i, od in enumerate(inputKeys["bondlength"], start=1):
                atom_indices = get_atom_indices(inputKeys["fragment_indices"], od["bond_definition"])
                atoms: List[Atom] = [complexMolecule[atom_indices[0]], complexMolecule[atom_indices[1]]]
                label = f"bondlength{i}_{'-'.join([f'{str(atom_index)}{atom.symbol}' for atom_index, atom in zip(atom_indices, atoms)])}"
                outputData[label] = atoms[0].distance_to(atoms[1]) - od["original_value"]

        if inputKeys["angle"]:
            for i, od in enumerate(inputKeys["angle"], start=1):
                atom_indices = get_atom_indices(inputKeys["fragment_indices"], od["angle_definition"])
                atoms: List[Atom] = [complexMolecule[atom_indices[0]], complexMolecule[atom_indices[1]], complexMolecule[atom_indices[2]]]
                label = f"angle{i}_{'-'.join([f'{str(atom_index)}{atom.symbol}' for atom_index, atom in zip(atom_indices, atoms)])}"
                outputData["angle"] = atoms[0].angle(atoms[1], atoms[2]) - od["original_value"]

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
