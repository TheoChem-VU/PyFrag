import pathlib as pl
import re
from typing import Any, Callable, Dict, List, Sequence, Tuple, Union

from pyfrag.errors import PyFragCoordFileError, PyFragSectionInputError
from pyfrag.read_input.pyfrag_settings import (
    Angle,
    BondLength,
    Dihedral,
    OrbitalEnergy,
    OrbitalEnergyWithIrrep,
    Overlap,
    OverlapWithIrrep,
    Population,
    PopulationWithIrrep,
    PyFragSection,
)

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


def _read_bondangle_line(line: str) -> Tuple[int, int, float]:
    """Reads the line containing the "angle" keyword. The correct format for the line is:

    angle atom1 atom2 [angle] (optional)

    """
    line_content: List[str] = _check_line_length(line, "angle", (3, 4))

    # Two atoms without angle
    if len(line_content) == 3:
        atom1, atom2 = line_content[1:]
        return int(atom1), int(atom2), 0.0

    # Two atoms with angle
    atom1, atom2, angle = line_content[1:]
    return int(atom1), int(atom2), float(angle)


def _read_dihedral_angle(line: str) -> Tuple[int, int, int, float]:
    """Reads the line containing the "dihedral" keyword. The correct format for the line is:

    dihedral atom1 atom2 atom3 [dihedral_angle] (optional)

    """
    line_content: List[str] = _check_line_length(line, "dihedral", (4, 5))

    # Three atoms without angle
    if len(line_content) == 4:
        atom1, atom2, atom3 = line_content[1:]
        return int(atom1), int(atom2), int(atom3), 0.0

    # Three atoms with angle
    atom1, atom2, atom3, dihedral_angle = line_content[1:]
    return int(atom1), int(atom2), int(atom3), float(dihedral_angle)


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

    or:

    fragment 1 2 5-8
    fragment -1 (this will select all atoms except the ones of the other fragment)

    """
    line_content: List[str] = _check_line_length(line, "fragment", (2, 100_000), strict_limit=False)
    _, unprocessed_indices = line_content[0], line_content[1:]

    if any(index == "0" for index in unprocessed_indices):
        raise PyFragSectionInputError("fragment is not valid. The atom indices should start with 1.", "fragment")

    if len(unprocessed_indices) == 1 and unprocessed_indices[0] == "-1":
        return [-1]

    indices: List[int] = []
    for index in unprocessed_indices:
        if "-" in index:
            start, end = index.split("-")
            indices.extend(list(range(int(start), int(end) + 1)))
        else:
            indices.append(int(index))

    return indices


def _read_coord_file_line(line: str) -> str:
    """Reads the line containing the "coord_file" keyword. Multiple coordinate files may be specified by having multiple lines. The correct format for the line is:

    coord_file FILENAME1
    coord_file FILENAME2

    """
    line_content: List[str] = _check_line_length(line, "coord_file", (2, 2))
    _, coord_files = line_content

    return coord_files


read_functions: Dict[str, Callable[[str], Any]] = {
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


def extract_pyfrag_section(pyfrag_section: str, input_file_dir: pl.Path) -> PyFragSection:
    """Extracts extra specifications from the PyFrag input file.

    This function takes the PyFrag input file as an argument and extracts extra specifications such as orbitalenergy, overlap, population of a certain fragment and MO.
    The function returns a dictionary containing the extracted specifications.

    Args:
        pyfrag_section (str): The PyFrag section of the input file. It must be a string with newline characters.

    Returns:
        dict[str, Any]: A dictionary containing the extracted specifications.

    Example of the returned dictionary:
    {
        'bondlength': [(1, 2, 1.0)],
        'angle': [(1, 2, 120.0)],
        'dihedral': [(1, 2, 3, 180.0)],
        'fragment': [[1, 2, 3, 4], [5, 6, 7, 8]]
        'strain': [-0.5, 0.5],
        'overlap': [('frag1', 'HOMO', 'frag2', 'LUMO')]
        'coord_file': file.xyz
    }

    """
    pyfrag_lines = pyfrag_section.split("\n")
    input_keys: Dict[str, list[Any]] = {}

    for line in pyfrag_lines:
        line = line.lower()

        for key, message in deprecated_keys_and_message_mapping.items():
            if line.startswith(key):
                raise PyFragCoordFileError(message)

        # Property keys such as bondlength, angle, dihedral, overlap, population, orbitalenergy, vdd
        for key, func in read_functions.items():
            if line.startswith(key):  # This checks for comments such as #, !, :
                input_keys[key] = [] if key not in input_keys else input_keys[key]
                input_keys[key].append(func(line))

    # Convert the input_keys dictionary to a PyFragSection model
    return PyFragSection(
        bondlengths=[BondLength(atom1=bl[0], atom2=bl[1], bond_length=bl[2]) for bl in input_keys.get("bondlength", [])],
        angles=[Angle(atom1=a[0], atom2=a[1], angle=a[2]) for a in input_keys.get("angle", [])],
        dihedrals=[Dihedral(atom1=d[0], atom2=d[1], atom3=d[2], dihedral_angle=d[3]) for d in input_keys.get("dihedral", [])],
        overlaps=[
            Overlap(frag1=o[0], MO1=o[1], frag2=o[2], MO2=o[3]) if len(o) == 4 else OverlapWithIrrep(irrep1=o[0], frag1=o[1], index1=o[2], irrep2=o[3], frag2=o[4], index2=o[5])
            for o in input_keys.get("overlap", [])
        ],
        populations=[Population(frag1=p[0], MO1=p[1]) if len(p) == 2 else PopulationWithIrrep(irrep1=p[0], frag1=p[1], index1=p[2]) for p in input_keys.get("population", [])],
        orbitalenergies=[OrbitalEnergy(frag1=oe[0], MO1=oe[1]) if len(oe) == 2 else OrbitalEnergyWithIrrep(irrep1=oe[0], frag1=oe[1], index1=oe[2]) for oe in input_keys.get("orbitalenergy", [])],
        vdd_indices=[item for sublist in input_keys.get("vdd", []) for item in sublist],  # Flatten the list of lists
        strain_values=[s for s in input_keys.get("strain", [])],
        fragment_indices=[f for f in input_keys.get("fragment", [])],
        coordfile=[input_file_dir / c for c in input_keys.get("coordfile", [])],
    )
