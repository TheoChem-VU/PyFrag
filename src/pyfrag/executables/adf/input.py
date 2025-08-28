import logging
import os
import pathlib as pl
import re
from typing import Any, Dict, List, Sequence, Tuple, TypeVar, Union

# Since we use amspython (from the AMS pacakge) as python environment, we cannot install the pyfrag package.
# However, when performing tests with pytest, we install the pyfrag package and then the "." needs to be added.
# Therefore, we include this try ... except block.
try:
    from .adf_parser import extract_sections
    from .errors import PyFragCoordFileError, PyFragSectionInputError
    from .mol_handling import find_files_with_wildcard_option
except ImportError:
    from adf_parser import extract_sections
    from errors import PyFragCoordFileError, PyFragSectionInputError
    from mol_handling import find_files_with_wildcard_option

logger = logging.getLogger("Input Reader-ADF")


# For Python 3.7, which is used by AMS2021 and older, we cannot use the TypedDict implementation since it's not supported (only in Python 3.8 and later)
try:
    from typing import TypedDict  # type: ignore # Nasty import but necessary for compatibility
except ImportError:
    # Minimal fallback for Python <3.8 without typing-extensions
    class TypedDict(dict):
        pass


def expandvars_backslash(path: pl.Path) -> pl.Path:
    """Short function to expand environment variables such as $SLURM_SUBMIT_DIR"""
    expanded_pathstring = re.sub(r"(?<!\\)\$[A-Za-z_][A-Za-z0-9_]*", "", os.path.expandvars(path))
    return pl.Path(expanded_pathstring)


# =============================================================================
# TypedDicts for input keys
# =============================================================================


class OrbitalDescription(TypedDict):
    type: str  # The type of the orbital energy (e.g. "INDEX", "HOMO", "LUMO")
    frag: str  # The fragment specifier (e.g. "frag1", "frag2")
    irrep: Union[str, None]  # The irreducible representation if applicable (only when the format "[irrep] [frag] [index]" is used)
    index: Union[str, None]  # The index of the orbital if applicable (only when the format "[irrep] [frag] [index]" is used)


class BondLength(TypedDict):
    atom_indices: List[int]  # The indices of the atoms involved in the bond
    original_value: float  # The original value of the bond length if applicable


class Angle(TypedDict):
    atom_indices: List[int]  # The indices of the atoms involved in the angle
    original_value: float  # The original value of the angle if applicable


class Dihedral(TypedDict):
    atom_indices: List[int]  # The indices of the atoms involved in the dihedral
    original_value: float  # The original value of the dihedral if applicable


class InputKeys(TypedDict):
    restart_dir_name: Union[str, None]  # The name of the restart directory
    job_name: str  # The name of the output file
    log_level: int  # The log level for the program
    coordFile: List[pl.Path]  # The coordinate file(s) to be considered for reading in the trajectory
    fragment_energies: Dict[str, float]  # The equilibrium energy for each fragment
    fragment_indices: Dict[str, List[int]]  # collection of fragment indices, e.g. {"frag1": [1, 2], "frag2": [3, 4]}
    bondlength: List[BondLength]  # Atom indices (and optionally the original value) for which bond lengths are printed
    angle: List[Angle]  # Atom indices (and optionally the original value) for which angles are printed
    dihedral: List[Dihedral]  # Atom indices (and optionally the original value) for which dihedrals are printed
    VDD: Sequence[int]  # Atom numbers for which VDD charges are printed
    hirshfeld: Sequence[str]  # Fragment specifiers (e.g. "frag1", "frag2") for which Hirshfeld charges are printed
    irrepOI: List[Dict[str, str]]  # Irreps for which Orbital Interactions (OI) energies are printed
    population: List[OrbitalDescription]  # Orbital populations for each fragment
    overlap: List[Tuple[OrbitalDescription, OrbitalDescription]]  # Overlap between two fragment orbitals
    orbitalenergy: List[OrbitalDescription]  # Orbital energies for each fragment
    adfinputfile: Union[str, None]  # ADF input file for the general settings
    old_adfinputfile: Union[str, None]  # ADF input file for the general settings (old version)
    fragment1_extra: Union[str, None]  # ADF input file for the first fragment (extra settings)
    fragment2_extra: Union[str, None]  # ADF input file for the second fragment (extra settings)
    complex_extra: Union[str, None]  # ADF input file for the complex (extra settings)


# =============================================================================
# Reading PyFrag Section functions
# =============================================================================

deprecated_keys_and_message_mapping: Dict[str, str] = {
    "ircpath": "'ircpath' is deprecated. Use 'coordfile' instead.",
}


def _check_line_length(line: str, input_key: str, limits: Sequence[int], strict_limit: bool = True) -> List[str]:
    """Checks if the line has the correct length for reading the specified input key.

    This function checks if the line has the correct length for reading the specified input key. If the line does not have the correct length, an error is raised.
    The function returns a list containing the values of the line.

    Args:
        line (str): The line containing the keyword and the values.
        input_key (str): The keyword to be read.
        limits (Sequence[int]): The limits of the length of the line such as (3, 4) for bondlength and angle.
        strict_limit (bool, optional): If True, the length of the line must be exactly the same as the limits. If False, the length of the line must be within the limits. Defaults to True.

    Raises:
        PyFragSectionInputError: If the line does not have the correct length.

    Returns:
        List[str]: A list containing the values of the line.

    """
    line_content: list[str] = re.split(r"\s*[#!;:]\s*", line.strip())[0].split()

    if strict_limit and len(line_content) not in limits:
        raise PyFragSectionInputError(f"Length of the {input_key} is not correct. Make sure to specify the correct format", input_key)
    elif not strict_limit and len(line_content) < limits[0] or len(line_content) > limits[1]:
        raise PyFragSectionInputError(f"Length of the {input_key} is not correct. Make sure to specify the correct format", input_key)
    return line_content


def _read_bondlength_line(line: str) -> Tuple[int, int, float]:
    """Reads the line containing the "bondlength" keyword. The correct format for the line is:

    bondlength atom1 atom2 [bondlength] (optional)

    """
    line_content: List[str] = _check_line_length(line, "bondlength", (3, 4))

    # Bondlength has not been specified
    if len(line_content) == 3:
        atom1, atom2 = line_content[1:]
        return int(atom1), int(atom2), 0.0

    # Bondlength has been specified
    atom1, atom2, bondlength = line_content[1:]
    return int(atom1), int(atom2), float(bondlength)


def _read_bondangle_line(line: str) -> Tuple[int, int, int, float]:
    """Reads the line containing the "angle" keyword. The correct format for the line is:

    angle atom1 atom2 atom3 [angle] (optional)

    """
    line_content: List[str] = _check_line_length(line, "angle", (4, 5))

    # Two atoms without angle
    if len(line_content) == 4:
        atom1, atom2, atom3 = line_content[1:]
        return (int(atom1), int(atom2), int(atom3), 0.0)

    # Two atoms with angle
    atom1, atom2, atom3, angle = line_content[1:]
    return int(atom1), int(atom2), int(atom3), float(angle)


def _read_dihedral_angle(line: str) -> Tuple[int, int, int, int, float]:
    """Reads the line containing the "dihedral" keyword. The correct format for the line is:

    dihedral atom1 atom2 atom3 atom4 [dihedral_angle] (optional)

    """
    line_content: List[str] = _check_line_length(line, "dihedral", (5, 6))

    # Three atoms without angle
    if len(line_content) == 5:
        atom1, atom2, atom3, atom4 = line_content[1:]
        return (int(atom1), int(atom2), int(atom3), int(atom4), 0.0)

    # Three atoms with angle
    atom1, atom2, atom3, atom4, dihedral_angle = line_content[1:]
    return int(atom1), int(atom2), int(atom3), int(atom4), float(dihedral_angle)


def _read_overlap_line(line: str) -> Union[Tuple[str, str, str, str, str, str], Tuple[str, str, str, str]]:
    """Reads the line containing the "overlap" keyword. The correct formats are:

    overlap frag1 HOMO frag2 LUMO
    overlap frag1 HOMO-1 frag2 LUMO+3
    overlap S frag1 5 AA frag2 4

    """
    line_content: List[str] = _check_line_length(line, "overlap", (5, 7))

    # Two fragments, two MOs (HOMO/LUMO kind), no irreps
    if len(line_content) == 5:
        frag1, MO1, frag2, MO2 = line_content[1:5]

        # Check if the fragments are strings and are named frag1 and frag2
        assert all(["frag" in frag for frag in [frag1, frag2]])
        return str(frag1), str(MO1), str(frag2), str(MO2)

    irrep1, frag1, index1, irrep2, frag2, index2 = line_content[1:]
    return str(irrep1), str(frag1), str(index1), str(irrep2), str(frag2), str(index2)


def _read_population_line(line: str) -> Union[Tuple[str, str], Tuple[str, str, str]]:
    """Reads the line containing the "population" keyword. The correct formats are:

    population frag1 HOMO
    population frag2 HOMO-1
    population AA frag2 5

    """
    line_content: List[str] = _check_line_length(line, "population", (3, 4))

    # Two fragments, two MOs (HOMO/LUMO kind), no irreps
    if len(line_content) == 3:
        (
            frag1,
            MO1,
        ) = line_content[1:3]
        return str(frag1), str(MO1)

    irrep1, frag1, index1 = line_content[1:]
    return str(irrep1), str(frag1), str(index1)


def _read_orbitalenergy_line(line: str) -> Union[Tuple[str, str], Tuple[str, str, str]]:
    """Reads the line containing the "orbitalenergy" keyword. The correct formats are:

    orbitalenergy frag1 HOMO
    orbitalenergy frag1 HOMO-2
    orbitalenergy AA frag2 5

    """
    line_content: List[str] = _check_line_length(line, "orbitalenergy", (3, 4))

    # Two fragments, two MOs (HOMO/LUMO kind), no irreps
    if len(line_content) == 3:
        (
            frag1,
            MO1,
        ) = line_content[1:3]
        return str(frag1), str(MO1)

    irrep1, frag1, index1 = line_content[1:]
    return str(irrep1), str(frag1), str(index1)


def _read_vdd_line(line: str) -> List[int]:
    """Reads the line containing the "vdd" keyword. The correct formats are:

    vdd 1 2 3 4
    vdd 3 6 8

    """
    line_content: list[str] = re.split(r"\s*[#!;:]\s*", line.strip())[0].split()
    try:
        all([int(atom_index) for atom_index in line_content[1:]])
    except ValueError:
        raise PyFragSectionInputError("Make sure to specify the vdd charges with spaces in between the indices", "vdd")

    return [int(atom_index) for atom_index in line_content[1:]]


def _read_irrep_line(line: str) -> List[str]:
    """Reads the line containing the "irrep" keyword. The correct formats are:

    irrepOI AA
    irrepOI E1:1

    """
    line_content: List[str] = _check_line_length(line, "orbitalenergy", (2, 2))

    irrep = line_content[1]
    return [irrep]


def _read_strain_line(line: str) -> float:
    """Reads the line containing the "strain" keyword. The correct formats are:

    strain  0.5
    strain -0.5

    """
    line_content: List[str] = _check_line_length(line, "strain", (2, 2))
    _, value = line_content

    try:
        value = float(value)
    except ValueError:
        raise PyFragSectionInputError("Make sure to specify the strain value correctly with a float number", "strain")

    return value


def _read_fragment_indices_line(line: str) -> List[int]:
    """Reads the line containing the "fragment" keyword which specifies the atom indices of each fragment. The correct formats are:

    fragment 1 2 4
    fragment 3 6 8

    or:

    fragment 1-4
    fragment 5-8

    """
    line_content: List[str] = _check_line_length(line, "fragment", (2, 100_000), strict_limit=False)
    _, unprocessed_indices = line_content[0], line_content[1:]

    if any(index == "0" for index in unprocessed_indices):
        raise PyFragSectionInputError("fragment is not valid. The atom indices should start with 1.", "fragment")

    indices: List[int] = []
    for index in unprocessed_indices:
        if "-" in index:
            start, end = index.split("-")
            indices.extend(list(range(int(start), int(end) + 1)))
        else:
            indices.append(int(index))

    return indices


def _read_coord_file_line(line: str) -> str:
    """Reads the line containing the "coord_file" keyword. Multiple coordinate files may be specified by giving a (space-separated) list. The correct format for the line is:

    coord_file FILENAME1 FILENAME2

    """
    line_content: List[str] = _check_line_length(line, "coord_file", (2, 100_000), strict_limit=False)
    _, coord_files = line_content

    return coord_files


# Define the mapping of line parsers
read_functions = {
    "bondlength": _read_bondlength_line,
    "angle": _read_bondangle_line,
    "dihedral": _read_dihedral_angle,
    "overlap": _read_overlap_line,
    "population": _read_population_line,
    "orbitalenergy": _read_orbitalenergy_line,
    "vdd": _read_vdd_line,
    "irrep": _read_irrep_line,
    "strain": _read_strain_line,
    "fragment": _read_fragment_indices_line,
    "coordfile": _read_coord_file_line,
}


def extract_pyfrag_section(pyfrag_section: str) -> Dict[str, Any]:
    """Extracts extra specifications from the PyFrag input file.

    This function takes the PyFrag input file as an argument and extracts extra specifications such as orbitalenergy, overlap, population of a certain fragment and MO.
    The function returns a dictionary containing the extracted specifications.

    Args:
        pyfrag_section (str): The PyFrag section of the input file. It must be a string with newline characters.
        input_file_dir (pl.Path): The directory containing the input file.

    Returns:
        Dict[str, Any]: A dictionary containing the extracted specifications.
    """
    pyfrag_lines = pyfrag_section.split("\n")
    input_keys: Dict[str, Any] = {}

    for line in pyfrag_lines:
        # Skip empty lines and comments
        if not line.strip() or line.strip().startswith("#"):
            continue

        line_lower = line.lower()

        # Replace all commas with spaces since commas are not allowed in the input
        line = line.replace(",", " ")

        # Handle simple key-value pairs
        if " " in line_lower:
            parts = line_lower.split()
            key = parts[0]

            # Handle specific keys with custom parsing
            if key in read_functions:
                input_keys[key] = input_keys.get(key, [])
                input_keys[key].append(read_functions[key](line))
            elif key in ["name", "job_name", "jobname"]:
                input_keys["job_name"] = parts[1]
            elif key == "log_level":
                input_keys["log_level"] = parts[1].upper()
            elif key.startswith("frag") and key.endswith("_indices"):
                # Handle frag1_indices, frag2_indices, etc.
                input_keys[key] = _read_fragment_indices_line(line)
            elif key.startswith("frag") and key.endswith("_energy"):
                # Handle frag1_energy, frag2_energy, etc.
                input_keys[key] = float(parts[1])
            else:
                # Default handling for other keys
                input_keys[key] = " ".join(parts[1:]) if len(parts) > 1 else ""

    return input_keys


# =============================================================================
# Main function to process user input
# =============================================================================

TList = TypeVar("TList", bound=List[Union[int, str, float]])


def _flatten_list(input_list: TList) -> TList:
    """Flatten nested lists into a single list."""
    flat_list: TList = []  # type: ignore  # This type definition is not correct since we make a new list insread of updating the input_list
    for item in input_list:
        if isinstance(item, list):
            flat_list.extend(_flatten_list(item))
        else:
            flat_list.append(item)
    return flat_list


def process_user_input(input_file_path: str) -> InputKeys:
    """
    Process user input from file and return InputKeys dictionary.

    Args:
        input_file_path (str): Path to the input file (like symmetry_job.in)

    Returns:
        InputKeys: Dictionary containing all parsed input parameters
    """
    # Extract sections from the input file
    input_file_sections = extract_sections(input_file_path)

    # Set default values
    inputKeys: InputKeys = InputKeys(
        restart_dir_name=None,
        job_name="PyFragJob",
        log_level=logging.INFO,
        fragment_energies={},
        fragment_indices={},
        coordFile=[],
        VDD=[],
        hirshfeld=[],
        bondlength=[],
        angle=[],
        dihedral=[],
        irrepOI=[],
        population=[],
        overlap=[],
        orbitalenergy=[],
        adfinputfile=None,
        old_adfinputfile=None,
        fragment1_extra=None,
        fragment2_extra=None,
        complex_extra=None,
    )

    # Process the PyFrag section if it exists
    if "pyfrag" in input_file_sections:
        pyfrag_data = extract_pyfrag_section(input_file_sections["pyfrag"])

        # Process the extracted data
        for key, val in pyfrag_data.items():
            key_lower = key.lower()

            # Handle overlap
            if key_lower == "overlap":
                inputValue = []
                for term in val:
                    if len(term) == 4:
                        inputValue.append((OrbitalDescription(type=term[1], frag=term[0], irrep=None, index=None), OrbitalDescription(type=term[3], frag=term[2], irrep=None, index=None)))
                    else:
                        inputValue.append((OrbitalDescription(type="INDEX", frag=term[1], irrep=term[0], index=term[2]), OrbitalDescription(type="INDEX", frag=term[4], irrep=term[3], index=term[5])))
                inputKeys["overlap"] = inputValue

            # Handle population
            elif key_lower == "population":
                inputValue = []
                for term in val:
                    if len(term) == 2:
                        inputValue.append(OrbitalDescription(type=term[1], frag=term[0], irrep=None, index=None))
                    else:
                        inputValue.append(OrbitalDescription(type="INDEX", frag=term[1], irrep=term[0], index=term[2]))
                inputKeys["population"] = inputValue

            # Handle orbitalenergy
            elif key_lower == "orbitalenergy":
                inputValue = []
                for term in val:
                    if len(term) == 2:
                        inputValue.append(OrbitalDescription(type=term[1], frag=term[0], irrep=None, index=None))
                    else:
                        inputValue.append(OrbitalDescription(type="INDEX", frag=term[1], irrep=term[0], index=term[2]))
                inputKeys["orbitalenergy"] = inputValue

            # Handle coordinate files
            elif key_lower in ["coordfile", "ircpath", "irct21", "lt"]:
                if key_lower != "coordfile":
                    logger.warning(f"Input option '{key_lower}' is deprecated. Please use 'coordfile' instead.")

                if isinstance(val, list):
                    inputKeys["coordFile"] = [coord_file for coord_file in val]
                else:
                    inputKeys["coordFile"] = [val]

                # First handle environment variables such as "$SLURM_SUBMIT_DIR"
                inputKeys["coordFile"] = [expandvars_backslash(path) for path in inputKeys["coordFile"]]

                # Then, expand the wildcard option
                inputKeys["coordFile"] = find_files_with_wildcard_option(inputKeys["coordFile"])

                # Check if all coordinate files exist.
                if not all(path.exists() for path in inputKeys["coordFile"]):
                    raise PyFragCoordFileError("One or more coordinate files do not exist. Please check the paths.", inputKeys["coordFile"])

            # Handle fragment indices
            elif key_lower == "fragment" or key_lower == "fragment_indices":
                for frag_index, indices in enumerate(val, start=1):
                    frag_key = f"frag{frag_index}"
                    inputKeys["fragment_indices"][frag_key] = indices

            # Handle specific fragment indices (frag1_indices, frag2_indices, etc.)
            elif re.match(r"frag\d+_indices", key_lower):
                match = re.match(r"frag(\d+)_indices", key_lower)
                if match:
                    frag_index = int(match.group(1))
                    inputKeys["fragment_indices"][f"frag{frag_index}"] = val

            # Handle fragment energies
            elif key_lower == "strain":
                inputKeys["fragment_energies"] = {f"frag{i}": energy for i, energy in enumerate(val, start=1)}

            elif re.match(r"frag\d+_energy", key_lower):
                match = re.match(r"frag(\d+)_energy", key_lower)
                if match:
                    frag_index = int(match.group(1))
                    inputKeys["fragment_energies"][f"frag{frag_index}"] = float(val)

            # Handle VDD
            elif key_lower == "vdd":
                inputKeys["VDD"] = _flatten_list(val)

            # Handle hirshfeld
            elif key_lower == "hirshfeld":
                inputKeys["hirshfeld"] = _flatten_list([frag for frag in val if frag.startswith("frag")])

            # Handle bondlength
            elif key_lower == "bondlength":
                inputValue = []
                for term in val:
                    if len(term) == 2:
                        inputValue.append(BondLength(atom_indices=[term[0], term[1]], original_value=0.0))
                    else:
                        inputValue.append(BondLength(atom_indices=[term[0], term[1]], original_value=term[2]))
                inputKeys["bondlength"] = inputValue

            # Handle angle
            elif key_lower == "angle":
                inputValue = []
                for term in val:
                    if len(term) == 3:
                        inputValue.append(Angle(atom_indices=[term[0], term[1], term[2]], original_value=0.0))
                    else:
                        inputValue.append(Angle(atom_indices=[term[0], term[1], term[2]], original_value=term[3]))
                inputKeys["angle"] = inputValue

            # Handle dihedral angles
            elif key_lower == "dihedral":
                inputValue = []
                for term in val:
                    if len(term) == 4:
                        inputValue.append(Dihedral(atom_indices=[term[0], term[1], term[2], term[3]], original_value=0.0))
                    else:
                        inputValue.append(Dihedral(atom_indices=[term[0], term[1], term[2], term[3]], original_value=term[4]))
                inputKeys["dihedral"] = inputValue

            # Handle ADF input files
            elif key_lower in ["adfinputfile", "old_adfinputfile", "fragment1_extra", "fragment2_extra", "complex_extra"]:
                inputKeys[key_lower] = val

            # Handle job name
            elif key_lower in ["job_name", "name", "jobname"]:
                inputKeys["job_name"] = val

            # Handle log level
            elif key_lower == "log_level":
                log_mapping = {"DEBUG": logging.DEBUG, "INFO": logging.INFO, "WARNING": logging.WARNING, "ERROR": logging.ERROR, "CRITICAL": logging.CRITICAL}
                inputKeys["log_level"] = log_mapping.get(val.upper(), logging.INFO)

            # Handle other keys
            else:
                if isinstance(val, list):
                    inputKeys[key_lower] = val
                else:
                    inputKeys[key_lower] = [val]

    # Copy other sections to InputKeys
    for section_name, section_content in input_file_sections.items():
        if section_name.lower() in ["fragment1_extra", "fragment2_extra", "complex_extra"]:
            inputKeys[section_name.lower()] = section_content
        elif section_name.lower() == "ams":
            inputKeys["adfinputfile"] = section_content
        elif section_name.lower() == "adf":
            inputKeys["old_adfinputfile"] = section_content

    return inputKeys
