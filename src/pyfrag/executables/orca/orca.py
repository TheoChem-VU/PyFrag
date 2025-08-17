import math
import os
import pathlib as pl
import re
import sys
import tempfile
from typing import Any, Dict, List, Optional, TextIO, Tuple, Union


def read_energy_t21(t21: str, unit: str = "kcal/mol") -> float:
    """
    Reads the total energy from a T21 file and returns it in kcal/mol.

    Args:
        t21 (str): Path to the T21 file.
        unit (str): Energy unit (default: "kcal/mol").

    Returns:
        float: Total energy in kcal/mol.
    """
    with open(t21, "r") as input_file:
        gout = input_file.read()
    # Match any line with 'energy' and 'Eh', extract the number before Eh
    match = re.search(r"energy.*?:\s*([-+]?\d*\.\d+|\d+)\s*Eh", gout, re.IGNORECASE)
    if not match:
        raise ValueError(f"Could not find energy value in {t21}")
    try:
        energy = float(match.group(1))
        return energy * 627.51
    except Exception as e:
        raise ValueError(f"Could not parse energy from line: {match.group(0)}") from e


def bondlength_xyz(xyz: List[List[float]], atom1: Union[int, str], atom2: Union[int, str]) -> Optional[float]:
    """
    Calculates the bond length between two atoms in a given xyz matrix.

    Args:
        xyz (List[List[float]]): Coordinates matrix.
        atom1 (Union[int, str]): Index of the first atom (1-based).
        atom2 (Union[int, str]): Index of the second atom (1-based).

    Returns:
        Optional[float]: Bond length or None if input is invalid.
    """
    try:
        bond = [int(atom1), int(atom2)]
    except Exception:
        sys.stderr.write(f"Values for printing bondlength are not integers!\nGot: {atom1}, {atom2}\n")
        return None

    try:
        bondlength = math.sqrt((xyz[bond[0] - 1][0] - xyz[bond[1] - 1][0]) ** 2 + (xyz[bond[0] - 1][1] - xyz[bond[1] - 1][1]) ** 2 + (xyz[bond[0] - 1][2] - xyz[bond[1] - 1][2]) ** 2)
        return bondlength
    except IndexError:
        sys.stderr.write(f"Atom indices out of range for bondlength calculation!\nGot: {bond[0]}, {bond[1]}\n")
        return None


def shelljoin(*args: str) -> str:
    """
    Joins path components and replaces backslashes with forward slashes.

    Args:
        *args (str): Path components.

    Returns:
        str: Joined path with forward slashes.
    """
    path = ""
    for a in args:
        path = os.path.join(path, a)
    path = re.sub(r"\\", "/", path)
    return path


def write_key(file: TextIO, value: Any, pform: str = r"%.4f", ljustwidth: int = 16) -> None:
    """
    Writes a value or list of values to a file, formatted for output.

    Args:
        file (TextIO): File object to write to.
        value (Any): Value or list of values to write.
        pform (str): Print format for floats.
        ljustwidth (int): Width for left-justification.
    """
    if not isinstance(value, list):
        value = [value]
    for val in value:
        if val is None:
            file.write(str.ljust(" None", ljustwidth))
        else:
            try:
                val_float = float(val)
                if val_float < 0:
                    val_str = str.ljust(pform % (val_float), ljustwidth)
                else:
                    val_str = " " + str.ljust(pform % (val_float), ljustwidth - 1)
                file.write(val_str)
            except Exception:
                file.write(str.ljust(str(val) + " ", ljustwidth - 1))


def print_comment_out(text: str, width: int = 80) -> None:
    """
    Prints a formatted comment block to stdout.

    Args:
        text (str): Comment text.
        width (int): Total width of the comment block.
    """
    comment = " ====================== COMMENT "
    comment = re.sub(r"COMMENT", text, comment)
    comment += "=" * (width - len(comment))
    comment = "\n" + comment + "\n"
    sys.stdout.write(comment)


def make_xyz_matrix(matrix: List[float], nr_of_atoms: int) -> List[List[List[float]]]:
    """
    Converts a flat list of coordinates into a 3D xyz matrix.

    Args:
        matrix (List[float]): Flat list of coordinates.
        nr_of_atoms (int): Number of atoms.

    Returns:
        List[List[List[float]]]: 3D xyz matrix.
    """
    bohr2angs = 1.88972612499
    xyz_matrix: List[List[List[float]]] = []
    for i in range(0, len(matrix) - nr_of_atoms, nr_of_atoms * 3):
        coor: List[List[float]] = []
        for j in range(0, nr_of_atoms * 3, 3):
            coor.append([matrix[i + j] / bohr2angs, matrix[i + j + 1] / bohr2angs, matrix[i + j + 2] / bohr2angs])
        xyz_matrix.append(coor)
    return xyz_matrix


def read_irc_xyz_out(main_dir: Union[str, pl.Path], output_file: str) -> None:
    """
    Reads coordinates from an IRC or LT output file and writes them to an .amv file.

    Args:
        main_dir (str): Main directory path.
        output_file (str): Output file name.
    """
    atomlist: Dict[int, str] = {1: "H", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 26: "Fe", 46: "Pd"}
    main_dir = pl.Path(main_dir).resolve()
    file_name = main_dir.joinpath(output_file)
    try:
        with open(file_name, "r") as input_file:
            gsscript = input_file.read()
    except Exception:
        sys.stderr.write("file " + str(file_name) + " could not be opened!\n\n")
        sys.exit(1)

    step = re.compile(r"# OF STEPS =" + r"\s*\d*", re.DOTALL).findall(gsscript)
    step = "".join(step)
    step = re.compile(r"# OF STEPS =", re.DOTALL).sub("", step)
    step = step.split()
    try:
        step = [int(x) for x in step]
    except Exception:
        raise ValueError("Could not parse step numbers from IRC output.")

    number = 0
    loop: List[int] = []
    for i in step:
        number += i
        loop.append(number)

    data: List[List[List[str]]] = []
    atomnumber: Optional[int] = None
    for i in loop:
        spec: List[List[str]] = []
        if re.compile(r"Input orientation:", re.IGNORECASE).search(gsscript):
            specs_found = re.compile(r"Input orientation:" + r".*?" + r"Distance matrix", re.IGNORECASE | re.DOTALL).findall(gsscript)
            if len(specs_found) > i:
                spec_str = specs_found[i]
                spec_str = re.compile(r"Distance matrix (angstroms):" + r".*", re.IGNORECASE).sub("", spec_str)
                spec_str = re.compile(r"Input orientation:" + r".*", re.IGNORECASE).sub("", spec_str)
                newspec: List[List[str]] = []
                lines = spec_str.splitlines()[5:-2]
                atomnumber = len(lines)
                for line in lines:
                    line_split = line.split()
                    newspec.append([atomlist.get(int(line_split[1]), "X"), line_split[3], line_split[4], line_split[5]])
                data.append(newspec)

    root, ext = os.path.splitext(output_file)
    ircfile = main_dir.joinpath(root + ".amv")
    with open(ircfile, "a") as f:
        for i in range(len(loop)):
            f.writelines("\n")
            if atomnumber is None:
                raise ValueError("Atom number could not be determined in IRC output.")
            for j in range(atomnumber):
                f.writelines("\n")
                test = data[i][j]
                f.writelines(["%s " % item for item in test])


def read_xyz_file(file_name: str) -> Tuple[List[str], List[List[List[float]]]]:
    """
    Reads coordinates from a (multi)-xyz file and returns atom list and xyz matrix.

    Args:
        file_name (str): Path to the xyz file.

    Returns:
        Tuple[List[str], List[List[List[float]]]]: Atom list and xyz matrix.
    """
    try:
        with open(file_name, "r") as f:
            f_lines = f.readlines()
    except Exception:
        sys.stderr.write("file " + file_name + " could not be opened!\n\n")
        sys.exit(1)

    conv_count = 0
    for l in f_lines:
        if re.compile(r"Geometry.*Converged", re.IGNORECASE).search(l):
            conv_count += 1

    xyz_matrix: List[List[List[float]]] = []
    xyz: List[List[float]] = []
    atom_list: List[str] = []
    xyz_line_found = False

    logfile = re.search(r"\.logfile", file_name)
    if conv_count < 2:
        logfile = False

    converged = True  # First geometry found in IRCs is the TS-geometry.

    for line in f_lines:
        if not converged and logfile:
            converged = bool(re.compile(r"Geometry.*Converged", re.IGNORECASE).search(line))
        elif not logfile:
            converged = True
        match = re.findall(r"\w{1,2}\b\s+.*\d+.*\s+.*\d+.*\s+.*\d.*\s", line)
        if match:
            temp_xyz = str.split(match[0])[1:4]
            try:
                for j in range(3):
                    temp_xyz[j] = float(temp_xyz[j])
                xyz.append(temp_xyz)
                if len(xyz_matrix) == 0 and converged:
                    atom_list.append(str.split(match[0])[0])
                xyz_line_found = True
            except Exception:
                xyz_line_found = False
        else:
            xyz_line_found = False
            if xyz and converged:
                xyz_matrix.append(xyz)
                xyz = []
                converged = False
            elif not converged:
                xyz = []

    if len(xyz) > 0:
        xyz_matrix.append(xyz)

    return atom_list, xyz_matrix


def new_adf(adfscript: str, dir: str, name: str, name_t21: str, atoms_block: str, extra: str, fragments: bool = False, title: str = "") -> None:
    """
    Creates and runs a new ADF script, writing input and output files.

    Args:
        adfscript (str): Basic ADF script.
        dir (str): Directory for files.
        name (str): Name for the file.
        name_t21 (str): Name for the T21 file.
        atoms_block (str): Atoms block for the script.
        extra (str): Extra input blocks.
        fragments (bool): Whether to use fragments.
        title (str): Title for the calculation.
    """
    new_adfscript = re.compile(r"end input", re.IGNORECASE).sub(atoms_block + "\n*" + "\n" + extra, adfscript)
    for _ in range(10):
        new_adfscript = re.sub(r"\n\n\n", "\n\n", new_adfscript)
    starting = """
 ====================== Starting ADF with inputscript: =========================

"""
    sys.stdout.write(starting)
    sys.stdout.write(new_adfscript)
    orca_input = name + ".inp"
    orca_output = shelljoin(dir, name + ".out")

    with open(orca_input, "w") as f:
        f.write("#orca block \n")
        f.write(new_adfscript)

    with open("orca_input_sh", "w") as f:
        f.write("orca " + orca_input + ">" + orca_output)

    os.system("sh orca_input_sh")
    os.remove("orca_input_sh")


def get_input_block(specs: str, match_string: str, match_string2: str, default: str = "") -> str:
    """
    Reads a block of input specifications from the input script.

    Args:
        specs (str): Input specification string.
        match_string (str): Start string to match.
        match_string2 (str): End string to match.
        default (str): Default value if nothing matches.

    Returns:
        str: The matched block or default.

    """
    spec = default
    if re.compile(match_string, re.IGNORECASE).search(specs):
        if match_string2 and re.compile(match_string, re.IGNORECASE).search(specs):
            spec = re.compile(match_string + r".*?" + match_string2, re.IGNORECASE | re.DOTALL).findall(specs)[0]
            spec = re.compile(match_string2 + r".*", re.IGNORECASE).sub("", spec)
            spec = re.compile(match_string + r".*", re.IGNORECASE).sub("", spec)
    return spec


def get_input_line(specs: str, match_string: str, default: Optional[Any] = None, dim: int = 0) -> Any:
    """
    Reads a line of input specifications from the input script.

    Args:
        specs (str): Input specification string.
        match_string (str): String to match.
        default (Optional[Any]): Default value if nothing matches.
        dim (int): Determines return type (0: str, 1: list, 2: 2D list).

    Returns:
        Any: Matched value or default.

    Examples:
        print bond 1 3 1.00
        print bond 1 2 2.00

        dim=2 -> [[1, 3, 1.00], [1, 2, 2.00]]
    """
    spec = default

    p = re.compile(match_string, re.IGNORECASE)
    spec = p.findall(specs)

    if not spec:
        return default

    # Check for "=" seperator because we only want the latter part (the "value", not the "key")
    # Also, remove white spaces because we don't want to create file such as 'pyfrag_dir/ frag1/...'
    for i in range(len(spec)):
        if "=" in spec[i]:
            spec[i] = spec[i].split("=")[1].strip()

        # Remove any whitespaces
        spec[i] = re.sub(r"\s+", " ", spec[i]).strip()

    # Unpack the (nested) list
    if dim == 0:
        if spec and isinstance(spec[0], list):
            spec = spec[0][0]
        else:
            spec = spec[0]
    elif dim == 1:
        # If the spec is a nested list, flatten it
        if isinstance(spec[0], list):
            spec = [item for sublist in spec for item in sublist]
        elif isinstance(spec[0], str):
            spec = spec[0].split()
    elif dim == 2:
        # If the spec is a nested list, flatten it
        if isinstance(spec[0], list):
            spec = [item for sublist in spec for item in sublist]
        elif isinstance(spec[0], str):
            spec = [line.split() for line in spec]

    return spec


class InputPrm:
    """
    Class representing an input parameter for PyFrag.
    """

    def __init__(self, name: str, specs: str = "add", match_string: str = "add", default: Optional[Any] = None, dim: int = 0, print_names: Union[bool, List[str]] = False, unit: str = "") -> None:
        self.name: str = name
        self.unit: str = unit
        self.args: Any = get_input_line(specs, match_string, default, dim)
        self.cv: Optional[Any] = None
        self.values: List[Any] = []
        self.columns: List[Any] = []

        if print_names:
            self.print_names: List[str] = print_names if isinstance(print_names, list) else [name]
        elif dim and self.args is not None:
            self.print_names = []
            for arg in self.args:
                temp_name = name
                for a in arg:
                    temp_name = temp_name + "_" + a
                self.print_names.append(temp_name)
        else:
            self.print_names = [name]

        self.pwidth: int = 16  # default columnwidth value

        for name_str in self.print_names:
            if len(name_str) > self.pwidth - 2:
                self.pwidth = len(name_str) + 2

        self.dim: int = dim

    def get(self) -> "InputPrm":
        return self

    def set(self, value: Any) -> None:
        self.cv = value

    def tell(self) -> None:
        """
        Print details of parameter.
        """
        print('Details of key with name "%s"' % self.name)
        print("nr of args:", self.args)
        print("names printed in data file", self.print_names)
        print("Dimension?:", self.dim)

    def __str__(self) -> str:
        return f"InputPrm(name={self.name}, args={self.args}, cv={self.cv}, values={self.values}, columns={self.columns}, print_names={self.print_names}, pwidth={self.pwidth}, dim={self.dim})"


class Parameters:
    """
    Class for managing a collection of InputPrm objects.
    """

    def __init__(self) -> None:
        self.dict: Dict[str, InputPrm] = {}
        self.seq: List[InputPrm] = []
        self.keys_list: List[str] = []

    def add(self, name: str, specs: str = "Default", match_string: str = "Default", default: Optional[Any] = None, dim: int = 0, print_names: Union[bool, List[str]] = False, unit: str = "") -> None:
        self.dict[name] = InputPrm(name, specs, match_string, default, dim, print_names, unit)
        self.seq.append(self.dict[name])
        self.keys_list.append(name)
        self.__dict__[name] = self.dict[name].get()

    def __getitem__(self, name: str) -> InputPrm:
        return self.dict[name]

    def __setitem__(self, name: str, value: Any) -> None:
        self.dict[name].set(value)

    def keys(self) -> List[str]:
        return list(self.dict.keys())

    def list_keys(self) -> List[str]:
        return self.keys_list

    def __str__(self) -> str:
        return "\n".join(str(param) for param in self.seq)


# -----------------------------------------------------#
#               M A I N  P R O G R A M                #
# -----------------------------------------------------#

header = """
 ===============================================================================
 *  Pyfrag 2019.02                                                             *
 *  Streamlining your reaction path analysis!                                  *
 *                                                                             *
 *  Author: Willem-Jan van Zeist                                               *
 *
 *                                                                             *
 *  Find the manual at http://www.few.vu.nl/~wolters/pyfrag/                   *
 *  For some examples on how to use PyFrag, see the examples directory.        *
 *                                                                             *
 *  E-mail for PyFrag: LP.Wolters@vu.nl                                        *
 ===============================================================================
"""
sys.stdout.write(header)
print_comment_out("PyFrag Is Initializing")

# input should be: pyfrag.py jobfile [tempdir]
if len(sys.argv) < 2:
    raise RuntimeError("No input file specified. Usage: pyfrag.py jobfile [tempdir]")
with open(sys.argv[1], "r") as input_file:
    adfscript = input_file.read()
input_file_dir = pl.Path(sys.argv[1]).parent

# separate input specifications and adf script
specs = re.compile(r"PYFRAG\s*INPUT.*END\s*PYFRAG\s*INPUT", re.DOTALL).findall(adfscript)
if not specs:
    specs = re.compile(r"INPUT_?\s*SPECS.*END INPUT_SPECS", re.DOTALL).findall(adfscript)
if specs:
    specs = specs[0]
else:
    sys.stderr.write("Pyfrag input statments not found! Exiting Pyfrag.\n")
    sys.exit(1)

adfscript = re.compile(r".*END\s*INPUT_?\s*SPECS|.*END\s*PYFRAG\s*INPUT", re.DOTALL).sub("", adfscript, 1)

# getting inline
inline = get_input_line(specs, r"inline\s(.*)", default=False)
if inline:
    try:
        with open(inline, "r") as inline_file:
            specs += inline_file.read()
    except Exception:
        sys.stderr.write("Inline file " + str(inline) + " could not be opened.\n")

if len(sys.argv) > 2:
    temp_dir = tempfile.mkdtemp(prefix="Pyfrag_TMPDIR_", dir=sys.argv[2])
else:
    temp_dir = tempfile.mkdtemp(prefix="Pyfrag_TMPDIR_")

# Reading whether fragment directory and/or datafile is specified
frag_dir = os.path.basename(get_input_line(specs, r"fragments?.*dir.(.*)", default="fragmentfiles"))
data_file = os.path.basename(get_input_line(specs, r"data.*file.(.*)", default="fragment_energies.txt"))

# check or make subdirectory
main_dir = input_file_dir  # this is the main directory for the energy_file
frag_dir = input_file_dir / frag_dir  # this is the dir for the fragment files

if not frag_dir.exists():
    frag_dir.mkdir(parents=True)

tape21_temp_dir = shelljoin(temp_dir, "TAPE21s")

try:
    os.mkdir(tape21_temp_dir)
except Exception:
    sys.stderr.write("Temporary TAPE21 dir is already present:" + tape21_temp_dir + "\n")

# Reading the user specs of inputfile
lt = get_input_line(specs, r"type\s*=\s*(LT)", default=None)
irc = get_input_line(specs, r"type\s*=\s*(IRC)", default=None)

if lt is None and irc is None:
    sys.stderr.write("The type of pyfrag calculation (IRC/LT/SP) has not been specified.\n")

# Reading all statements which indicate the printing of values --> printing parameters pp
pp = Parameters()

pp.add("Filename", specs, r"fa1?_name\s*=?(.*)", default="fa", unit=".out")
pp.Filename.pwidth = len(pp.Filename.args) + 7  # '_' and four numbers are added after the filename, increase column width (plus 2 default extra)
pp.add("Point", default=[], unit="Nr.")
pp.add("Bondlength", specs, r"print.*bond\s*=?\s*(\d\s*\d.*)", dim=2, unit="Angstrom")
pp.add("TotalIntEn", unit="kcal/mol")
pp.add("StrainFrag1", specs, r"print\s*strain\s*frag\s*1(.*)", unit="kcal/mol")
pp.add("StrainFrag2", specs, r"print\s*strain\s*frag\s*2(.*)", unit="kcal/mol")
if pp.StrainFrag1.args is not None and pp.StrainFrag2.args is not None:
    pp.add("StrainTotal", unit="kcal/mol")
    pp.add("EnergyTotal", unit="kcal/mol")

print_comment_out("Read print statements")

# Extra, if not present, they will be set as empty string and always passed on to the adfscript
extra_frag1 = get_input_block(specs, r"EXTRA\s*frag1", r"END\s*EXTRA\s*frag1")
extra_frag2 = get_input_block(specs, r"EXTRA\s*frag2", r"END\s*EXTRA\s*frag2")
extra_fa = get_input_block(specs, r"EXTRA\s*fa1?", r"END\s*EXTRA\s*fa")

# Reads specified names for fragments and fa-files
name_frag1 = get_input_line(specs, r"\s*frag1\s*=(.*)", default="frag1")
name_frag2 = get_input_line(specs, r"\s*frag2\s*=(.*)", default="frag2")

name_fa = get_input_line(specs, r"fa1?_name\s*=(.*)", default="fa")
name_fa_t21 = shelljoin(tape21_temp_dir, name_fa + ".t21")

# Defining the list of xyz_matrices for all inputs, if the frag_lists are not used (as in sp) they will be tested as false
matrix: List[List[List[float]]] = []
matrix_list: List[str] = []
matrix_list_frag1: List[str] = []
matrix_list_frag2: List[str] = []

# Read name for LT/IRC outputfile, list for the fragment atoms. One of the types should be specified.
# calls function to read all coordinates from the file
if not lt and not irc:
    raise ValueError("Neither the 'IRC' or 'LT' type has been specified in the 'type = []' line.")

output_file = get_input_line(specs, r"\s*output\s*file\s*1?\s*=?(.*)")
output_file_2 = get_input_line(specs, r"\s*output\s*file\s*2\s*=?(.*)")

if (output_file is None and output_file_2 is None) or not any(file is None for file in [output_file, output_file_2]):
    sys.stderr.write("Outputfile(s) for irc or lt are not specified\n")

LT_frag1 = get_input_block(specs, r"frag1\s*=", r"end frag1")
LT_frag1 = re.findall(r"\b\d+\b", LT_frag1)
LT_frag1 = [int(x) for x in LT_frag1]

LT_frag2 = get_input_block(specs, r"frag2\s*=", r"end frag2")
LT_frag2 = re.findall(r"\b\d+\b", LT_frag2)
LT_frag2 = [int(x) for x in LT_frag2]

# retrieving atom_list from input
atom_list2 = get_input_block(specs, r"frag1\s*=", r"end frag2")
atom_list2 = re.findall(r"\b\d+\b\s*\b[a-zA-Z]{1,2}\b", atom_list2)
for i in range(len(atom_list2)):
    atom_list2[i] = re.split(" ", atom_list2[i])
    atom_list2[i][0] = int(atom_list2[i][0])
atom_list2.sort()
for i in range(len(atom_list2)):
    atom_list2[i] = atom_list2[i][1]


atom_list = []
if any(ext in output_file for ext in ["xyz", "logfile", "amv"]):
    atom_list, matrix = read_xyz_file(os.path.join(main_dir, output_file))
    if output_file_2:
        atom_list, matrix2 = read_xyz_file(os.path.join(main_dir, output_file_2))
        matrix.reverse()
        matrix.extend(matrix2)

if any(ext in output_file for ext in ["log", "out"]):
    read_irc_xyz_out(main_dir, output_file)
    root, ext = os.path.splitext(output_file)
    atom_list, matrix = read_xyz_file(os.path.join(main_dir, root + ".amv"))
    if output_file_2:
        read_irc_xyz_out(main_dir, output_file_2)
        root, ext = os.path.splitext(output_file_2)
        atom_list, matrix2 = read_xyz_file(os.path.join(main_dir, root + ".amv"))
        matrix.reverse()
        matrix.extend(matrix2)

print("\n*** Input parameters ***")
print(pp)

number_of_loops = len(matrix)
alist_length = len(atom_list)
if len(atom_list2) < alist_length:
    alist_length = len(atom_list2)
for i in range(alist_length):
    if not re.compile(atom_list2[i], re.IGNORECASE).search(atom_list[i]):
        sys.stderr.write("\nyour range of atoms does not match that of the given output file!\n")
        sys.stderr.write("Your atoms:\n")
        for j in range(len(atom_list2)):
            sys.stderr.write(atom_list2[j] + "\n")
        sys.stderr.write("Atoms from outputfile:\n")
        for j in range(len(atom_list)):
            sys.stderr.write(atom_list[j] + "\n")
        sys.stderr.write("\nWill continue with the atoms with atomnumbers as in the inputfile!\n")

for xyz in matrix:
    fa_block, frag1_block, frag2_block = "", "", ""
    for j in range(len(xyz)):
        if j + 1 not in LT_frag1 and j + 1 not in LT_frag2:
            continue
        fa_block += atom_list[j] + "  "
        for k in range(len(xyz[j])):
            fa_block += str(xyz[j][k]) + "  "
        if j + 1 in LT_frag1:
            fa_block += " \n"
        if j + 1 in LT_frag2:
            fa_block += " \n"
    fa_block += "\n"
    matrix_list.append(fa_block)
    for nr in LT_frag1:
        frag1_block += atom_list[nr - 1] + "  "
        for k in range(len(xyz[nr - 1])):
            frag1_block += str(xyz[nr - 1][k]) + " "
        frag1_block += "\n"
    frag1_block += "\n"
    matrix_list_frag1.append(frag1_block)
    for nr in LT_frag2:
        frag2_block += atom_list[nr - 1] + "  "
        for k in range(len(xyz[nr - 1])):
            frag2_block += str(xyz[nr - 1][k]) + " "
        frag2_block += "\n"
    frag2_block += "\n"
    matrix_list_frag2.append(frag2_block)

# For both SP and IRC/LT: get the bondlengths, angles, etc.
if pp.Bondlength.args is not None:
    for i in range(len(pp.Bondlength.args)):
        bondlengths: List[Optional[float]] = []
        if len(pp.Bondlength.args[i]) > 2:
            diffb = float(pp.Bondlength.args[i][-1])
        else:
            diffb = 0
        for xyz in matrix:
            bondlength = bondlength_xyz(xyz, pp.Bondlength.args[i][0], pp.Bondlength.args[i][1])
            if bondlength is not None:
                bondlengths.append(bondlength - diffb)
            else:
                bondlengths.append(None)
        pp.Bondlength.args[i] = bondlengths

calc_points: List[int] = []

data_file_initialized = False

if not calc_points:
    for i in range(0, number_of_loops):
        calc_points.append(i)

pp.Point.args = []
for i in range(number_of_loops):
    pp.Point.args.append(i)
pp.Point.pwidth = len("Point") + 2
if pp.Point.pwidth < len(str(number_of_loops)) + 2:
    pp.Point.pwidth = len(str(number_of_loops)) + 2

print_comment_out("PyFrag Initilization done")
print_comment_out("Starting to run the ADF jobs")

# giant loop for both sp and irc/lt, for the latter it makes the fragment inputfiles
for i in calc_points:
    print_comment_out("Starting cycle: " + str(i))
    if matrix_list_frag1:
        new_name_frag1 = name_frag1 + "_" + (4 - len(str(i))) * "0" + str(i)
        new_name_frag1_t21 = shelljoin(frag_dir, new_name_frag1 + ".out")
        new_name_frag1_t21_fa = os.path.basename(new_name_frag1 + ".out")
        new_adf(adfscript, frag_dir, new_name_frag1, new_name_frag1_t21, matrix_list_frag1[i], extra_frag1, title="frag1")
    else:
        new_name_frag1_t21 = name_fa_t21
        new_name_frag1_t21_fa = name_fa_t21

    if matrix_list_frag2:
        new_name_frag2 = name_frag2 + "_" + (4 - len(str(i))) * "0" + str(i)
        new_name_frag2_t21 = shelljoin(frag_dir, new_name_frag2 + ".out")
        new_name_frag2_t21_fa = os.path.basename(new_name_frag2 + ".out")
        new_adf(adfscript, frag_dir, new_name_frag2, new_name_frag2_t21, matrix_list_frag2[i], extra_frag2, title="frag2")
    else:
        new_name_frag2_t21 = name_fa_t21
        new_name_frag2_t21_fa = name_fa_t21

    if pp.StrainFrag1.args is not None:
        pp.StrainFrag1.cv = read_energy_t21(new_name_frag1_t21) - float(pp.StrainFrag1.args)
    else:
        strain_frag1_new = False

    if pp.StrainFrag2.args is not None:
        pp.StrainFrag2.cv = read_energy_t21(new_name_frag2_t21) - float(pp.StrainFrag2.args)
    else:
        strain_frag2_new = False

    new_name_fa = name_fa + "_" + (4 - len(str(i))) * "0" + str(i)
    new_name_fa_t21 = shelljoin(frag_dir, new_name_fa + ".out")

    new_adf(adfscript, frag_dir, new_name_fa, new_name_fa_t21, matrix_list[i], extra_fa, title="fragment_analysis")

    # Initializing data file if applicable
    if not data_file_initialized:
        energy_file = open(os.path.join(main_dir, data_file), "w")
        # Print headers of columns
        for key in pp.list_keys():
            if pp[key].args is not None:
                for name in pp[key].print_names:
                    energy_file.write(str.ljust(name, pp[key].pwidth))
        energy_file.write("\n")
        for key in pp.list_keys():
            if pp[key].args is not None:
                for name in pp[key].print_names:
                    energy_file.write(str.ljust(pp[key].unit, pp[key].pwidth))
        energy_file.write("\n")
        energy_file.close()

    data_file_initialized = True

    pp.Filename.cv = new_name_fa

    if pp.Bondlength.args is not None:
        pp.Bondlength.cv = []
        for bond in range(len(pp.Bondlength.args)):
            pp.Bondlength.cv.append(pp.Bondlength.args[bond][i])

    pp.Point.cv = pp.Point.args[i]

    print_comment_out("PyFrag is gathering data to print")

    # pp = get_print_values(new_name_fa_t21, pp)

    if pp.StrainFrag1.args is not None and pp.StrainFrag2.args is not None:
        pp.StrainTotal.cv = pp.StrainFrag1.cv + pp.StrainFrag2.cv
        pp.EnergyTotal.cv = read_energy_t21(new_name_fa_t21) - float(pp.StrainFrag1.args) - float(pp.StrainFrag2.args)
        pp.TotalIntEn.cv = pp.EnergyTotal.cv - pp.StrainTotal.cv
    energy_file = open(os.path.join(main_dir, data_file), "a")
    for key in pp.list_keys():
        if pp[key].args is None:
            continue
        if key == "Point":
            write_key(energy_file, pp[key].cv, pform="%.0f", ljustwidth=pp[key].pwidth)
            continue
        write_key(energy_file, pp[key].cv, ljustwidth=pp[key].pwidth)
        pp[key].values.append(pp[key].cv)
    energy_file.write("\n")

    energy_file.close()

    print_comment_out("PyFrag finished writing to data file")

print_comment_out("PyFrag finished. Have a nice day!")
