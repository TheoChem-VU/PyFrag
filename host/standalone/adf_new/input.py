import logging
import pathlib as pl
import re
from typing import Dict, List, Sequence, Tuple, Union

try:
    from typing import TypedDict
except ImportError:
    from typing_extensions import TypedDict

logger = logging.getLogger("PyFrag Parser")


class OrbitalDefition(TypedDict, total=True):
    type: str  # The type of the orbital energy (e.g. "INDEX", "HOMO", "LUMO")
    frag: str  # The fragment specifier (e.g. "frag1", "frag2")
    irrep: Union[str, None]  # The irreducible representation if applicable (only when the format "[irrep] [frag] [index]" is used)
    index: Union[str, None]  # The index of the orbital if applicable (only when the format "[irrep] [frag] [index]" is used)


class BondLength(TypedDict, total=True):
    bond_definition: List[int]  # The indices of the atoms involved in the bond
    original_value: float  # The original value of the bond length if applicable


class Angle(TypedDict, total=True):
    angle_definition: List[int]  # The indices of the atoms involved in the angle
    original_value: float  # The original value of the angle if applicable


class InputKeys(TypedDict, total=True):
    jobstate: Union[str, None]  # The name of the restart directory
    job_name: str  # The name of the output file
    log_level: int  # The log level for the program
    fragment_energies: Dict[str, float]  # The equilibrium energy for each fragment
    coordFile: List[pl.Path]  # The coordinate file(s) to be considered for reading in the trajectory
    fragment_indices: Dict[str, List[int]]  # collection of fragment indices, e.g. {"frag1": [1, 2], "frag2": [3, 4]}
    VDD: Sequence[int]  # Atom numbers for which VDD charges are printed
    hirshfeld: Sequence[str]  # Fragment specifiers (e.g. "frag1", "frag2") for which Hirshfeld charges are printed
    bondlength: List[BondLength]  # Atom indices (and optionally the original value) for which bond lengths are printed
    angle: List[Angle]  # Atom indices (and optionally the original value) for which angles are printed
    irrepOI: List[Dict[str, str]]  # Irreps for which Orbital Interactions (OI) energies are printed
    population: List[OrbitalDefition]  # Orbital populations for each fragment
    overlap: List[Tuple[OrbitalDefition, OrbitalDefition]]  # Overlap between two fragment orbitals
    orbitalenergy: List[OrbitalDefition]  # Orbital energies for each fragment
    adfinputfile: Union[str, None]  # ADF input file for the general settings
    old_adfinputfile: Union[str, None]  # ADF input file for the general settings (old version)
    fragment1_extra: Union[str, None]  # ADF input file for the first fragment (extra settings)
    fragment2_extra: Union[str, None]  # ADF input file for the second fragment (extra settings)
    complex_extra: Union[str, None]  # ADF input file for the complex (extra settings)


def process_user_input(parser) -> InputKeys:
    """
    Process user input and update the inputKeys dictionary.
    """
    # Set default values
    inputKeys: InputKeys = InputKeys(
        jobstate=None,
        job_name="PyFragJob",
        log_level=logging.DEBUG,
        fragment_energies={},
        fragment_indices={},
        coordFile=[],
        VDD=[],
        hirshfeld=[],
        bondlength=[],
        angle=[],
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

    for key, val in vars(parser.parse_args()).items():
        if val is not None:
            key = str(key).lower()
            inputValue = []
            # Possible options for overlap: "frag1 HOMO" = "[frag] [type]" or "AA frag1 2" = "[irrep] [frag] [index]" with [type] = "INDEX"
            if key == "overlap":
                for term in val:
                    if len(term) == 4:
                        # inputValue.append(({"type": term[1], "frag": term[0]}, {"type": term[3], "frag": term[2]}))
                        inputValue.append((OrbitalDefition(type=term[1], frag=term[0], irrep=None, index=None), OrbitalDefition(type=term[3], frag=term[2], irrep=None, index=None)))
                    else:
                        # inputValue.append(({"type": "INDEX", "frag": term[1], "irrep": term[0], "index": term[2]}, {"type": "INDEX", "frag": term[4], "irrep": term[3], "index": term[5]}))
                        inputValue.append((OrbitalDefition(type="INDEX", frag=term[1], irrep=term[0], index=term[2]), OrbitalDefition(type="INDEX", frag=term[4], irrep=term[3], index=term[5])))
                inputKeys[key] = inputValue

            elif key == "population":
                for term in val:
                    if len(term) == 2:
                        # inputValue.append({"type": term[1], "frag": term[0]})
                        inputValue.append(OrbitalDefition(type=term[1], frag=term[0], irrep=None, index=None))
                    else:
                        # inputValue.append({"type": "INDEX", "frag": term[1], "irrep": term[0], "index": term[2]})
                        inputValue.append(OrbitalDefition(type="INDEX", frag=term[1], irrep=term[0], index=term[2]))
                inputKeys[key] = inputValue

            elif key == "orbitalenergy":
                for term in val:
                    if len(term) == 2:
                        # inputValue.append({"type": term[1], "frag": term[0]})
                        inputValue.append(OrbitalDefition(type=term[1], frag=term[0], irrep=None, index=None))
                    else:
                        # inputValue.append({"type": "INDEX", "frag": term[1], "irrep": term[0], "index": term[2]})
                        inputValue.append(OrbitalDefition(type="INDEX", frag=term[1], irrep=term[0], index=term[2]))
                inputKeys[key] = inputValue

            elif key == "irrep OI" or key == "irrepOI" or key == "irrep" or key == "irrep_oi" or key == "irrep oi":
                logger.warning("The key 'irrep OI' is deprecated. The irreps are automatically detected and printed. Please remove this key from your input file.")
                # inputKeys["irrepOI"] = [{"irrep": term[0]} for term in val]

            elif key in ["ircpath", "irct21", "lt", "coordfile"]:
                if key != "coordfile":
                    logger.warning(f"Input option '{key}' is deprecated. Please use 'coordfile' instead.")

                # If there are multiple coordinate files, append them to the list
                inputKeys["coordFile"] = [pl.Path(term[0]) for term in val]

            elif key == "fragment" or key == "fragment_indices":
                for frag_index, indices in enumerate(val, start=1):
                    frag_key = f"frag{frag_index}"
                    inputKeys["fragment_indices"][frag_key] = [(i) for i in indices]

            # Match "frag[index]indices" pattern for the fragment indices
            elif re.match(r"frag\d+_indices", key) is not None:
                match = re.match(r"frag(\d+)_indices", key)
                if match:
                    frag_index = int(match.group(1))
                    inputKeys["fragment_indices"][f"frag{frag_index}"] = val[0]

            elif key == "strain":
                inputKeys["fragment_energies"] = {"frag" + str(i): a[0] for i, a in enumerate(val, start=1) if len(a) == 1}

            # Match "frag[index]_energy" pattern for the fragment energies
            elif re.match(r"frag\d+_energy", key) is not None:
                match = re.match(r"frag(\d+)_energy", key)
                if match:
                    frag_index = int(match.group(1))  # Extract the fragment index
                    inputKeys["fragment_energies"][f"frag{frag_index}"] = float(val[0])

            elif key == "VDD":
                inputKeys[key] = val[0]

            elif key == "hirshfeld":
                inputKeys[key] = [fragment_specifier for fragment_specifier in val[0] if fragment_specifier.startswith("frag")]

            elif key == "bondlength":
                for term in val:
                    if len(term) == 2:
                        inputValue.append((BondLength(bond_definition=[int(term[0]), int(term[1])], original_value=0.0)))
                    else:
                        inputValue.append((BondLength(bond_definition=[int(term[0]), int(term[1])], original_value=float(term[2]))))
                inputKeys[key] = inputValue

            elif key == "angle":
                for term in val:
                    if len(term) == 3:
                        inputValue.append((Angle(angle_definition=[int(term[0]), int(term[1]), int(term[2])], original_value=0.0)))
                    else:
                        inputValue.append((Angle(angle_definition=[int(term[0]), int(term[1]), int(term[2])], original_value=float(term[3]))))
                inputKeys[key] = inputValue

            elif key == "adfinputfile" or key == "old_adfinputfile" or key == "fragment1_extra" or key == "fragment2_extra" or key == "complex_extra":
                inputKeys[key] = val[0]

            elif key == "job_name" or key == "name" or key == "jobname":
                inputKeys["job_name"] = val

            elif key == "log_level":
                value = val[0].upper()
                log_mapping = {"DEBUG": logging.DEBUG, "INFO": logging.INFO, "WARNING": logging.WARNING, "ERROR": logging.ERROR, "CRITICAL": logging.CRITICAL}
                log_level = log_mapping.get(value, logging.INFO)

                if log_level is None:
                    raise ValueError(f"Invalid log level: {value}. Valid options are 'debug', 'info', 'warning', 'error', 'critical'.")
            else:
                inputKeys[key] = [term for term in val]

    return inputKeys
