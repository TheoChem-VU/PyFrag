import logging
import pathlib as pl
import re
import shutil
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple, Union

import constants as const
from errors import FragmentOptimizationError, PyFragSortComplexMoleculeError
from input import InputKeys
from mol_handling import create_pyfrag_trajectory_from_coord_file, update_fragment_indices
from result_classes import get_pyfrag_results
from scm.plams import AMSJob, AMSResults, Molecule, Settings, load_all

logger = logging.getLogger("PyFragModules-ADF")


# =====================================================================
# Handling of input files and settings
# =====================================================================


def settings_from_ams_block(ams_block: str) -> Settings:
    """
    Returns a Settings object from an AMS block

    Adapted from the AMSJob.from_inputfile method to read in the settings from a file and return a Settings object
    Reason for adapting was because the AMSJob.from_inputfile method does remove the ams block from the settings object
    """

    # Temporarily make an input file with solely the ams input block due to PLAMS2021 compatibility (it does not have the "from_input" method)
    with pl.Path("ams_input.tmp").open("w") as f:
        f.write(ams_block)
        f.flush()
        pre_calc_job: AMSJob = AMSJob.from_inputfile("ams_input.tmp")
    pl.Path("ams_input.tmp").unlink()  # Remove the temporary file after reading it

    # This does not include the "System" block; it does include the "Task" and "Engine" block
    settings = pre_calc_job.settings

    # Add "Task" block if it is not present
    if "Task" not in settings.input.ams:
        settings.input.ams.Task = "SinglePoint"

    # Work around for the fact that the ams block is not included in the settings object
    # The above code apparently does not allow for overriding the molecule
    if pre_calc_job.molecule is not None:
        molecule: Molecule = pre_calc_job.molecule[""]  # type: ignore
        charge = molecule.properties.charge
        settings.input.ams.System.Charge = charge

    return settings


# =====================================================================
# Handling of restarting the PyFrag Calculations
# =====================================================================


def handle_restart(foldername: Union[str, pl.Path], overwrite: bool = True) -> Optional[str]:
    """
    Handles restarting by checking if the folder exists, renaming it to a backup folder with a ".res" suffix,
    and returning the backup folder name. Supports both string and pathlib.Path inputs.

    Input:
       foldername (Union[str, pl.Path]): The foldername to restart from.
       overwrite (bool): Whether to overwrite the backup folder if it exists.

    Returns:
         Optional[str]: The name of the backup folder if it exists, otherwise None.
    """
    foldername = pl.Path(foldername).resolve()  # Ensure foldername is a pathlib.Path object
    restart_backup = None

    if foldername.is_dir():
        if any(foldername.iterdir()):  # Check if the folder is not empty
            restart_backup = foldername.with_suffix(".res")

            if overwrite and restart_backup.exists():
                shutil.rmtree(restart_backup)
                foldername.rename(restart_backup)
                logger.info(f"RESTART: Overwriting and moving {foldername.name} to {restart_backup.name}, and restarting from it")
            else:
                n = 1
                while restart_backup.exists():
                    n += 1
                    restart_backup = foldername.with_suffix(f".res{n}")
                foldername.rename(restart_backup)
                logger.info(f"RESTART: Moving {foldername.name} to {restart_backup.name} and restarting from it")
    else:
        logger.info("RESTART: The folder specified for restart does not exist, starting from scratch")

    return str(restart_backup) if restart_backup else None


# =====================================================================
# Handling of optimising the fragments if no strain energy is specified for the fragments
# =====================================================================


def update_fragment_strain_energies(input_keys: "InputKeys", fragment_jobs: Sequence[AMSJob]) -> "InputKeys":
    """Update the strain energies of the fragments in the input keys dictionary."""
    for i, job in enumerate(fragment_jobs, start=1):
        if job.results.ok():
            input_keys["fragment_energies"][f"frag{i}"] = job.results.get_energy(unit="kcal/mol")
        else:
            raise FragmentOptimizationError(f"Fragment {i} optimization failed. Please check the log and out file.")
    return input_keys


def optimize_fragments(frag1_mol: Molecule, frag2_mol: Molecule, frag1Settings: Settings, frag2Settings: Settings) -> List[AMSJob]:
    """Optimizes the fragments if requested by the user. Returns a list of AMSResults with a length of two."""
    job_names = ["frag1_opt", "frag2_opt"]
    opt_jobs = []
    for frag_mol, frag_settings in zip([frag1_mol, frag2_mol], [frag1Settings, frag2Settings]):
        frag_settings.input.ams.Task = "GeometryOptimization"
        frag_settings.input.ams.GeometryOptimization.Convergence.Gradients = "1e-4"
        job = AMSJob(molecule=frag_mol, settings=frag_settings, name=job_names.pop(0))
        opt_jobs.append(job)
    [job.run() for job in opt_jobs]
    [clean_up_job_folder(job) for job in opt_jobs]
    return opt_jobs


# =====================================================================
# Helper function for keeping the user-defined atom order of the fragments the same as in the complex molecule
# This is critical for ensuring the correct mapping of atoms between fragments and the original molecule, and thus extracting results from atom and geometric properties (e.g., charges, bond lengths, angles, etc.)
# =====================================================================


def sort_molecule_by_indices(molecule: Molecule, indices: List[int]) -> Molecule:
    """Sorts the molecule by the given indices which is used to map the atoms from the additional of the two fragments to the original molecule"""

    mol_copy = molecule.copy()
    n_atoms = len(mol_copy.atoms)
    indices_copy = indices.copy()

    # Convert to 0-based indexing if indices are 1-based
    zero_based_indices = [idx - 1 if idx > 0 else idx for idx in indices_copy]

    if any(index >= n_atoms or index < 0 for index in zero_based_indices):
        raise PyFragSortComplexMoleculeError(f"Indices {indices_copy} are out of bounds for the molecule with {n_atoms} atoms.")

    # Create a mapping from original position to desired position
    index_map = {original_idx: new_idx for new_idx, original_idx in enumerate(zero_based_indices)}

    # Sort the atoms based on the desired order specified by the fragment atom indices definitions (i.e., frag1_indices and frag2_indices keys)
    mol_copy.atoms = sorted(mol_copy.atoms, key=lambda atom: index_map.get(mol_copy.atoms.index(atom), float("inf")))

    return mol_copy


# =====================================================================
# Writing the results (see PyFrag section) in table format to a file
# =====================================================================

# Dealing with the units ----------------------------------------------

key_to_print_unit_mapping: Dict[FrozenSet[str], str] = {
    frozenset({"#IRC"}): "---",
    frozenset({"EnergyTotal", "Int", "Elstat", "Pauli", "OI", "Disp"}): "kcal/mol",
    frozenset({"StrainTotal", "frag1Strain", "frag2Strain"}): "kcal/mol",
    frozenset({"bondlength"}): "Angstrom",
    frozenset({"angle", "dihedral"}): "degrees",
    frozenset({"VDD"}): "millielectrons",
    frozenset({"hirshfeld"}): "electrons",
    frozenset({"orbitalenergy"}): "eV",
    frozenset({"population"}): "electrons",
    frozenset({"overlap"}): "---",
}


def get_unit_belonging_to_header_key(header: str) -> str:
    """
    Get the unit belonging to a specific header key.
    For example, "EnergyTotal" would return "kcal/mol".
    """
    for key_set, unit in key_to_print_unit_mapping.items():
        if any(header == key or header.startswith(key) for key in key_set):
            return unit
    return ""


# Dealing with ordering the header (keys) naturally----------------------------


def find_coordinates_axis(headers: Sequence[str]) -> Union[str, None]:
    """
    Find the coordinate axis on which the EDA is plotted on. It scans for the first header that starts with "bondlength", or "angle" if it exists.
    """
    for header in headers:
        if header.startswith("bondlength"):
            return header
        elif header.startswith("angle"):
            return header

    return None


def natural_sort_key(key: str) -> list:
    """
    Generate a natural sort key for strings containing numbers.
    For example, "bondlength10_xxx" will be sorted after "bondlength2_xxx", so the order will be: bondlength1_xxx, bondlength2_xxx, ..., bondlength10_xxx, bondlength11_xxx, ...
    The key may contain other numbers after the underscore which should be ignored.
    """
    if "_" in key:
        key = key.split("_")[0]
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", key)]


def write_table(data_rows: List[Dict[str, Union[str, float]]], output_file_name: str):
    """
    Write the results table to a file.
    There are certain keys that are always printed first so we need to order them first.
    Thereafter, the remaining keys are sorted naturally (from 1 to n).

    The table is formatted as:
    - Header row
    - Unit row
    - Data row1
    - Data row2
    - ...

    Or visually,

    |   #IRC   | EnergyTotal |    Int     |  Elstat   |  Pauli   |   OI   |  Disp   | StrainTotal | frag1Strain | frag2Strain | other keys |
    |----------|-------------|------------|-----------|----------|--------|---------|-------------|-------------|-------------|------------|
    |          |  kcal/mol   |  kcal/mol  | kcal/mol  | kcal/mol |kcal/mol|kcal/mol |  kcal/mol   |  kcal/mol   |  kcal/mol   |   [unit]   |
    |    1     |     ...     |     ...    |    ...    |   ...    |  ...   |   ...   |     ...     |    ...      |    ...      |    ...     |
    |    2     |     ...     |     ...    |    ...    |   ...    |  ...   |   ...   |     ...     |    ...      |    ...      |    ...     |
    """
    logger.info("Writing PyFrag results | EDA/ASM terms in kcal/mol | Orbital energies in eV | Bondlengths in Angstrom | (Dihedral) angles in degrees | VDD charges in millielectrons | Hirshfeld charges in electrons")
    standard_headers = ["#IRC", "EnergyTotal", "Int", "Elstat", "Pauli", "OI", "Disp", "StrainTotal", "frag1Strain", "frag2Strain"]

    coordinate_axis = find_coordinates_axis(list(data_rows[0].keys()))

    all_headers_sorted = sorted(data_rows[0], key=natural_sort_key)

    if coordinate_axis is not None:
        standard_headers.insert(1, coordinate_axis)

    selected_headers = [header for header in all_headers_sorted if header not in standard_headers]
    headers = standard_headers + selected_headers
    unit_row = [get_unit_belonging_to_header_key(header) for header in headers]

    # Calculate column widths dynamically based on the maximum length of headers and data
    column_widths = [max(len(header), max(len(str(data_row.get(header, ""))) for data_row in data_rows)) + 2 for header in headers]
    column_widths = [max(width, 9) for width in column_widths]  # Ensure minimum width of 9 to have short headers, but longer float data. This depends on the pform formatting in the write_row function

    with open(f"pyfrag_{output_file_name}.txt", "w") as output_file:
        # Write headers
        write_row(output_file, headers, ljustwidths=column_widths)
        write_row(output_file, unit_row, ljustwidths=column_widths)

        # Write data rows
        for data_row in data_rows:
            sorted_data_row = [data_row.get(header, "---") for header in headers]
            write_row(output_file, sorted_data_row, ljustwidths=column_widths)


def write_row(file, value, ljustwidths: Sequence[int], pform=r"%+7.4f"):
    # Write all data into a file with dynamic column widths
    for val, width in zip(value, ljustwidths):
        if val is None:
            file.write(str.ljust("  ---", width))
        else:
            if isinstance(val, float):
                file.write(str.ljust(pform % (val), width))
            else:
                file.write(str.ljust(str(val), width))
    file.write("\n")


def print_table(cellList, widthlist, bar):
    line = "-" * (sum(widthlist) + 4 * len(widthlist) + 6)
    if bar:
        print("\n", line)
    for i, entry in enumerate(cellList):
        print("  " + str(entry).ljust(widthlist[i]) + "  ", end=" ")
    print("")
    if bar:
        print(line)


# =====================================================================
# Cleaning up the calculation folders (complex.xxx and fragment[x].xxx folders) to remove unnecessary files (and save disk space)
# =====================================================================


def clean_up_job_folder(job: AMSJob):
    """
    Removes unnecessary files from the calculation folder
    See https://www.scm.com/doc/plams/components/results.html#scm.plams.core.results.Results._clean for more info
    """
    r: AMSResults = job.results
    if r.ok():
        mol = r.get_main_molecule()
        r._clean(["-", "$JN.err", "$JN.run", "CreateAtoms.out", "t12.rel"])
        for atom in set(mol.atoms):
            r._clean(["-", f"t21.*.{atom.symbol}"])
        job.pickle()  # this will update the .dill file which is used to restart the job and extract results when using plams


# =====================================================================
# Main function to run the PyFrag calculations
# =====================================================================


def pyfrag_driver(inputKeys: "InputKeys", frag1Settings: Settings, frag2Settings: Settings, complexSettings: Settings) -> Tuple[List[Dict[str, Union[str, float]]], "InputKeys"]:
    """
    Main function to run the PyFrag calculations. It will read in the input keys and settings, optimize the fragments if needed, and run the calculations for each point in the trajectory.
    It will return a dictionary with the updated input keys (such as the fragment indices and strain energies).
    """
    # main pyfrag driver used for fragment and complex calculation.
    # read coordinates from IRC or LT t21 file. Other choice is xyz file generated from other tools.
    load_all(inputKeys["restart_dir_name"]) if inputKeys["restart_dir_name"] is not None else None

    resultsTable = []

    molecule_trajectory: List[List[Molecule]] = create_pyfrag_trajectory_from_coord_file(coord_file=inputKeys["coordFile"], fragment_indices=list(inputKeys["fragment_indices"].values()))
    # update the fragment indices in the settings to match the indices in the trajectory

    updated_fragment_indices: List[List[int]] = update_fragment_indices(fragment_indices=list(inputKeys["fragment_indices"].values()))
    # Transform the list of lists into a dictionary with "frag1" and "frag2" as keys
    inputKeys["fragment_indices"] = {f"frag{i + 1}": frag_indices for i, frag_indices in enumerate(updated_fragment_indices)}

    length_of_trajectory: int = len(molecule_trajectory[0])

    # =====================================================================
    # Deciding if fragments need to be optimized.
    # If so, then add these optimized energies to the corresponding "strain" entries in the InputKeys which will be later used to calculate the total strain energy.
    # =====================================================================

    # Optimize fragments if the strain energy of both or one of the fragments is not given
    logger.info(msg="Checking if fragments need to be optimized")
    if len(inputKeys["fragment_energies"]) != 2:
        logger.info(msg="Optimizing fragments")
        frag1_mol, frag2_mol = molecule_trajectory[1][0], molecule_trajectory[2][0]
        optimized_frag_jobs = optimize_fragments(frag1_mol, frag2_mol, frag1Settings.copy(), frag2Settings.copy())  # copy settings to avoid changing the original settings
        logger.info("Updating strain energies of fragments")
        update_fragment_strain_energies(inputKeys, optimized_frag_jobs)
    else:
        logger.info(msg="Fragment strain energies are given, no need to optimize fragments")

    # =====================================================================
    # The main calculation loop over each point of the coordinate file
    # =====================================================================

    for path_index in range(1, len(molecule_trajectory[0]) + 1):
        logger.info(msg=f"Starting calculations for IRC point {path_index}/{length_of_trajectory}")
        # Iterate over all the first entries complex, frag1, frag2) with the second entry given by path_index (the point in the trajectory)
        # The first entry is the complex, the second entry is frag1 and the third entry is frag2
        for system_name, molecule in zip(const.SYSTEM_NAMES, [trajectory[path_index - 1] for trajectory in molecule_trajectory]):
            logger.debug(msg=f"{system_name}: {molecule}")

        outputData = {}
        outputData["StrainTotal"] = 0

        fragment_settings = [frag1Settings.copy(), frag2Settings.copy()]
        frag_jobs: List[AMSJob] = []

        # =====================================================================
        # First consider the fragments, starting from frag1 and frag2
        # NOTE: Default fragment names are currently used which are the second and third elements in the const.SYSTEM_NAMES list
        # NOTE: This may change later (>PyFrag2025) if the fragment names can be given as input options. \
        # NOTE: However, this is not as straight-forward to implement since some keys depend on the fragment names. An example is the FragOccupations key in which the names of the fragments must be named explicitly.
        # =====================================================================

        for frag_index, base_fragment_name in enumerate(const.SYSTEM_NAMES[1:], start=1):
            # Goes from frag1 -> frag1.xxxx1
            fragment_name = f"{base_fragment_name}.{str(path_index).zfill(5)}"

            frag_mol: Molecule = molecule_trajectory[frag_index][path_index - 1]
            frag_job = AMSJob(molecule=frag_mol.copy(), settings=fragment_settings[frag_index - 1], name=fragment_name)
            frag_job.run()

            # check if the calculation is successful and log error message if not
            if not frag_job.results.ok():
                logger.critical(msg=f"Fragment calculation for {fragment_name} failed with error message: {frag_job.get_errormsg()}")

            outputData[base_fragment_name + "Strain"] = frag_job.results.get_energy(unit="kcal/mol") - inputKeys["fragment_energies"][base_fragment_name]
            outputData["StrainTotal"] += outputData[base_fragment_name + "Strain"]
            clean_up_job_folder(frag_job)
            frag_jobs.append(frag_job)

        # =====================================================================
        # Second, consider the complex, which is the first element in the const.SYSTEM_NAMES list
        # Prepare the fragment by placing the "adf.rkf=..." suffix after the atom coordinates
        # Also, link the fragment names to the fragment files with the tuple: (frag_job, "adf")
        # =====================================================================

        for frag_index, frag_job in enumerate(frag_jobs, start=1):
            if frag_job.molecule is not None:
                for at in frag_job.molecule:
                    at.properties.suffix = f"adf.f=frag{frag_index}"  # type: ignore  # properties is a Settings instance which does not have explicit type hints

        complexMolecule = frag_jobs[0].molecule.copy() + frag_jobs[1].molecule.copy()  # type: ignore  # The addition is always between two molecules
        # Resort complexMolecule to match the original molecule (the one that is used to create the trajectory)
        complexMolecule = sort_molecule_by_indices(complexMolecule, updated_fragment_indices[0] + updated_fragment_indices[1])

        complexSettings.input.adf.fragments.frag1 = (frag_jobs[0], "adf")
        complexSettings.input.adf.fragments.frag2 = (frag_jobs[1], "adf")
        jobComplex = AMSJob(molecule=complexMolecule, settings=complexSettings, name=f"{const.SYSTEM_NAMES[0]}.{str(path_index).zfill(5)}")
        logger.info(msg=f"Running complex {path_index}")
        jobComplex.run()

        if not jobComplex.results.ok():
            logger.critical(msg="Complex calculation has failed, please check your input settings")

        # =====================================================================
        # Finally, get the results of the complex calculation *if* it was successful
        # =====================================================================

        if jobComplex.ok():
            # collect all data and put it in a list
            outputData["#IRC"] = str(path_index)
            # collect all data that need to be printed
            pyfragResult = get_pyfrag_results(jobComplex, inputKeys)
            # convert multiple value into a dict. NOTE: here we take the original complex molecule as the above addition (frag1_mol + frag2_mol) reorders the atoms
            outputLine = pyfragResult.get_output_data(complexMolecule, outputData, inputKeys)
            # collect updated informaiton of each point calculation and print it on screen

            headerList = sorted(outputLine.keys())
            a = headerList.pop(headerList.index("#IRC"))
            b = headerList.pop(headerList.index("EnergyTotal"))
            headerList = [a, b] + headerList

            valuesList = [str(outputLine[i]) for i in headerList]
            widthlist = [max(len(str(valuesList[_])), len(str(headerList[_]))) for _ in range(len(valuesList))]
            print_table(headerList, widthlist, False)
            print_table(valuesList, widthlist, False)
            resultsTable.append(outputLine)
            clean_up_job_folder(jobComplex)
            logger.info(msg=f"IRC point {path_index}/{length_of_trajectory} finished")
        else:
            logger.critical(msg=f"Complex calculation for {path_index} failed with error message: {jobComplex.get_errormsg()}")

    if len(resultsTable) == 0:
        raise RuntimeError("Calculations for all points failed, please check your input settings")

    # Returns the results as one big table including the (updated) InputKeys. The main script (`adf.py`) uses the table to print it to a file.
    return resultsTable, inputKeys
