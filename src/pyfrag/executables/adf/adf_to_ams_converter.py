# Original author: SCM developers
# Modified by: Siebe Lekanne Deprez; PhD student of the TheoCheM group at the VU (Amsterdam)

import sys

from scm.input_parser import InputParser  # noqa: F401  # type: ignore  # This is part of the *pre-compiled* AMS python environment and cannot be installed via pip install
from scm.plams import AMSJob, Atom, Molecule, Settings


def main_converter(content: str, is_file: bool = False) -> Settings:
    """
    Module designed to convert an ADF inputfile (before 2019) to the new AMS inputfile format (>2019)
    Argument:
        file (str): path to the (old) ADF inputfile
    Returns:
        settings (plams Settings object): the converted AMS inputfile
    """
    print()
    print("# ==================================================")
    print("# Automatic conversion of ADF-2019 input to AMS-2020")
    print("# ==================================================")
    print()

    if is_file:
        with open(content, "r") as f:
            adf_input = "".join(f.readlines()).strip()
    else:
        adf_input = content

    # In case one calls "$ADFBIN/adf -n 1 << eor", as it sometimes done for create runs...
    adf_input = adf_input.lstrip("1")
    adf_input += "\n"

    if "inline " in adf_input.lower():
        sys.exit("ERROR: INLINE key found in the ADF input. Cannot automatically convert input files containing the INLINE key.")

    # Convert the text input to a settings object:
    input_parser = InputParser()
    adf_settings = input_parser.to_settings("pre_amsified_adf", adf_input)

    # Convert the settings object from the ADF input to a settings with the new AMS input:
    adf_to_ams = ADFToAMS(adf_settings)

    # Return the converted settings
    return adf_to_ams.return_converted_settings()


class ADFToAMS(object):
    def __init__(self, adf_settings):
        self.adf = adf_settings

        self.warnings = []
        self.notes = []

        self.notes.append('"TAPE21" is now called "adf.rkf" and is located in AMS results folder (by default "ams.results")')

        # Get the molecule. This will remove the 'Atoms' input block from self.adf
        self.molecule = self._get_molecule()

        # Do all adf to ams conversion:
        self.sett = Settings()
        self._convert_input()

        # let PLAMS generate the new AMS input:
        job = AMSJob(settings=self.sett, molecule=self.molecule)
        self.input = job.get_input()

    def is_new_input_valid(self):
        # run the parser to check that the generated input is valid:
        input_parser = InputParser()
        tmp = input_parser.to_settings("ams", self.input)
        if len(self.sett.input.adf) > 0 and len(tmp.adf) == 0:
            return False
        else:
            return True

    def return_converted_settings(self):
        if not self.is_new_input_valid():
            print("")
            print("---------------------------------------------------------------")
            print("WARNING! the generated input is not valid. See the error above!")
            print("---------------------------------------------------------------")
            print("")

        if self.warnings:
            print("# === WARNINGS ===")
            for w in self.warnings:
                print("# -", w)
            print("")

        if self.notes:
            print("# === NOTES ===")
            for n in self.notes:
                print("# -", n)

        return self.sett

    # ==============
    # Private stuff:
    # ==============

    def _get_molecule(self):
        """
        Extract and return the molecule from the 'atoms' block of an adf input.
        This will also strip the parsed sections from the adf.
        """

        def parse_atoms_line(line, unit):
            """
            Given a line from the 'atoms' block of an old ADF input, return the corresponding Atom object
            """
            atom_properties = Settings()
            at = line.split()

            # atom line can start with a number (e.g. '1. C 0 0 0') Get rid of the number:
            if at[0].replace(".", "").isdigit():
                at = at[1:]

            # atomic symbol can start with a number (e.g. '1.C 0 0 0'). Get rid of the number:
            if at[0].split(".")[0].isdigit():
                at[0] = at[0].split(".")[1]

            # it might be a ghost atom (e.g. Gh.C 0 0 0)
            if "gh." in at[0].lower():
                atom_properties.ghost = True
                at[0] = at[0].split(".", 1)[1]

            # atomic symbol can have a name (e.g. 'C.bla 0 0 0')
            if "." in at[0]:
                at[0], atom_properties.name = at[0].split(".", 1)

            # atom line might have some suffix info (e.g. 'C 0 0 0 f=my_fragment bla=blu')
            if len(at) > 4:
                suffix = " ".join(at[4:])
                suffix = suffix.replace("/", "|")
                suffix = suffix.replace("f=", "adf.f=")
                suffix = suffix.replace("F=", "adf.F=")
                suffix = suffix.replace("r=", "adf.r=")
                suffix = suffix.replace("R=", "adf.R=")
                atom_properties.suffix = suffix

            symbol = at[0]
            coords = [float(x) for x in at[1:4]]

            atom = Atom(symbol=symbol, coords=coords, unit=unit)
            atom.properties = atom_properties

            return atom

        molecule = Molecule()

        # 'ZMat' and 'Internal' coordinates conversion not supported:
        if "atoms" in self.adf:
            if "_h" in self.adf.atoms:
                header = self.adf.atoms._h.lower()
                if "zmat" in header or "z-matrix" in header or "z-mat" in header:
                    self.warnings.append('Cannot automatically convert the "Atoms" block in "z-matrix" format.')
                    del self.adf.atoms
                    return molecule
                if "internal" in header:
                    self.warnings.append('Cannot automatically convert the "Atoms" block in "internal" format.')
                    del self.adf.atoms
                    return molecule

        unit = "angstrom"

        if "units" in self.adf:
            self.warnings.append('The key "Units" was defined in the input. The atomic coordinates have been converted to the proper unit, but other geometrical data (if present) has not been converted.')
            if "length" in self.adf.units:
                unit = self.adf.units.length
            del self.adf.units

        if "atoms" in self.adf:
            for line in self.adf.atoms._1:
                molecule.add_atom(parse_atoms_line(line, unit))
            del self.adf.atoms

        if "create" in self.adf:
            symbol = self.adf.create.split()[0]
            self.adf.create = self.adf.create.replace("/atomicdata/", "/atomicdata/ADF/")
            molecule.add_atom(Atom(symbol=symbol, coords=[0.0, 0.0, 0.0]))

        return molecule

    def _convert_input(self):
        self.sett.input.ams.task = "SinglePoint"

        self._handle_spin_polarization_and_charge()
        self._handle_symmetry()
        self._handle_relativity()
        self._handle_basis()
        self._handle_properties()
        self._misc()
        self._handle_geometry_block()
        self._remove_keys_and_raise_warning()

        self.sett.input.adf += self.adf

    def _handle_spin_polarization_and_charge(self):
        spin_polarization = 0

        if "charge" in self.adf:
            val = self.adf["charge"].split()
            charge = val[0]
            if len(val) > 1:
                spin_polarization = val[1]

            self.sett.input.ams.System.Charge = charge
            del self.adf.charge

        if "unrestricted" in self.adf:
            self.sett.input.adf.Unrestricted = "Yes"

            self.sett.input.adf.SpinPolarization = spin_polarization
            del self.adf.unrestricted

    def _handle_symmetry(self):
        if "symmetry" not in self.adf or str.lower(self.adf.symmetry) == "auto":
            # self.sett.input.ams.System.Symmetrize = 'Yes'
            # self.sett.input.ams.Symmetry.SymmetrizeTolerance = '0.001'
            # self.notes.append('Unlike ADF2019, AMS does not symmetrize the structure by default. See "System -> Symmetrize" in the AMS driver manual.')
            # self.notes.append('The AMS default symmetrization tolerance is larger than the ADF2019 one. See "Symmetry -> SymmetrizeTolerance" in the AMS driver manual.')
            self.notes.append("Symmetrize is disabled for PyFrag calculations due to transformation errors")
        if "symmetry" in self.adf and "tol" in self.adf.symmetry:
            del self.adf.symmetry
            self.warnings.append('The "Symmetry" tolerance was not automatically converted. You should be use the new "SymmetryTolerance" key in ADF.')

    def _handle_relativity(self):
        self.notes.append("Relativity is standard scalar-ZORA now. See 'Relativity' in the AMS driver manual.")
        if "relativistic" in self.adf:
            rel_string = self.adf.relativistic.lower()
            self.sett.input.adf.Relativity.Level = "scalar"

            if "spin-orbit" in rel_string or "spinorbit" in rel_string:
                self.sett.input.adf.Relativity.Level = "Spin-Orbit"

            if "zora" in rel_string:
                self.sett.input.adf.Relativity.Formalism = "ZORA"
            elif "ra-x2c" in rel_string:
                self.sett.input.adf.Relativity.Formalism = "RA-X2C"
            elif "x2c" in rel_string:
                self.sett.input.adf.Relativity.Formalism = "X2C"
            else:
                self.sett.input.adf.Relativity.Formalism = "Pauli"

            if "sapa" in rel_string:
                self.sett.input.adf.Relativity.Potential = "SAPA"

            del self.adf.relativistic
        else:
            self.sett.input.adf.Relativity.Level = "scalar"
            self.notes.append("Scalar relativistic effects (ZORA) are included by default in the 2020 version of ADF.")

        if "noncollinear" in self.adf:
            self.sett.input.adf.Relativity.SpinOrbitMagnetization = "NonCollinear"
            del self.adf.noncollinear

        if "collinear" in self.adf:
            self.sett.input.adf.Relativity.SpinOrbitMagnetization = "Collinear"
            del self.adf.collinear

        if "souexact" in self.adf:
            self.sett.input.adf.Relativity.souexact = "Yes"
            del self.adf.souexact

        if "soux" in self.adf:
            self.sett.input.adf.Relativity.SpinOrbitMagnetization = "CollinearX"
            del self.adf.soux

        if "souy" in self.adf:
            self.sett.input.adf.Relativity.SpinOrbitMagnetization = "CollinearY"
            del self.adf.souy

    def _handle_properties(self):
        if "analyticalfreq" in self.adf:
            self.sett.input.ams.Properties.NormalModes = "Yes"
            del self.adf.analyticalfreq

        if "gradient" in self.adf:
            self.sett.input.ams.Properties.Gradients = "Yes"
            del self.adf.gradient

        if "vcd" in self.adf:
            self.sett.input.ams.Properties.VCD = "Yes"
            del self.adf.vcd

        if "vroa" in self.adf:
            self.sett.input.ams.Properties.VROA = "Yes"
            del self.adf.vroa

    def _remove_keys_and_raise_warning(self):
        removed_message = 'The key "{}" does not exist anymore and it has been removed from the input. '

        to_be_removed = [
            ("Atomprops", 'See "System -> Atoms": "mass" and "adf.nuclear_charge" in the AMS driver manual.'),
            ("Constraints", 'See the "Constraints" key in the AMS driver manual.'),
            ("Restraint", 'See the "Restraints" key in the AMS driver manual.'),
            ("Geovar", "Check your input!"),
            ("TSRC", 'See the "TransitionStateSearch -> ReactionCoordinate" key in the AMS driver manual.'),
            ("Thermo", 'See the "Thermo" key in the AMS driver manual.'),
            ("QMMM", "To perform QMMM calculations you should use the new Hybrid engine in AMS."),
            ("Vibron", 'See "Task -> VibrationalAnalysis" in the AMS driver manual.'),
            ("RamanRange", 'See "Properties -> Raman" and "Raman -> FreqRange" in the AMS driver manual.'),
            ("ScanFreq", 'See "Properties -> NormalModes" and "NormalModes -> ReScanFreqRange" in the AMS driver manual.'),
            ("IRCStart", 'See "Task -> IRC" and "IRC" in the AMS driver manual.'),
            ("Noupdfac", "The Hessian update is now part of AMS."),
            ("BondOrder", 'See "Properties -> BondOrders" in the AMS driver manual.'),
        ]

        for key in [
            "database",
            "crdfilexyz",
            "crdfilemol",
            "cdatafile",
            "crdfilemopac",
            "grad_trf_btrf",
            "Hessdiag",
            "quild_nocoords_in_log",
            "sicoep",
            "lintermvxc",
            "nodrop",
            "prtiao",
            "readfcfile",
            "reducedmass",
            "sfguess",
            "solv",
            "testaor",
            "testfit",
            "testjob",
            "testpscharge",
            "trustsfguess",
            "userlegmn",
        ]:
            to_be_removed.append((key, ""))

        for key, warning in to_be_removed:
            if key in self.adf:
                del self.adf[key]
                self.warnings.append(removed_message.format(key) + warning)

    def _misc(self):
        if "aoresponse" in self.adf and "frequency" in self.adf.aoresponse:
            del self.adf.aoresponse.frequency
            self.warnings.append('The key "AOResponse -> Frequency" cannot be automatically converted, and has been removed from the input. You should new key "AOResponse -> Frequencies".')

        if "aoresponse" in self.adf and "freqrange" in self.adf.aoresponse:
            del self.adf.aoresponse.freqrange
            self.warnings.append('The key "AOResponse -> FreqRange" cannot be automatically converted, and has been removed from the input. You should new key "AOResponse -> Frequencies".')

        if "aoresponse" in self.adf and "efg" in self.adf.aoresponse:
            del self.adf.aoresponse.efg
            self.warnings.append('The key "AOResponse -> EFG" cannot be automatically converted, and has been removed from the input. You should use the new block "AOResponse -> EFG".')

        if "efield" in self.adf:
            self.sett.input.ams.System.ElectrostaticEmbedding.ElectricField = self.adf.efield + " [a.u.]"
            del self.adf.efield

        if "pointcharges" in self.adf:
            self.sett.input.ams.System.ElectrostaticEmbedding.MultipolePotential.Coordinates = self.adf.pointcharges
            del self.adf.pointcharges

        if "response" in self.adf and "frqbeg" in self.adf.response:
            del self.adf.response.frqbeg
            self.warnings.append('The key "Response -> Frqbeg" has been removed from the input. You should new key "Response -> Frequencies".')

        if "response" in self.adf and "frqend" in self.adf.response:
            del self.adf.response.frqend
            self.warnings.append('The key "Response -> Frqend" has been removed from the input. You should new key "Response -> Frequencies".')

        if "response" in self.adf and "nfreq" in self.adf.response:
            del self.adf.response.nfreq
            self.warnings.append('The key "Response -> Nfreq" has been removed from the input. You should new key "Response -> Frequencies".')

        if "fullscf" in self.adf:
            del self.adf.fullscf

        if "mblockbig" in self.adf:
            del self.adf.mblockbig

        if "extendpopan" in self.adf:
            del self.adf.extendpopan

        if "restart" in self.adf and "noexc" in self.adf.restart:
            del self.adf.restart.noexc

        if "restart" in self.adf and "nogeo" in self.adf.restart:
            del self.adf.restart.nogeo

        if "restart" in self.adf and "nohes" in self.adf.restart:
            del self.adf.restart.nohes

        if "restart" in self.adf and "file" in self.adf.restart:
            self.sett.input.ams.EngineRestart = self.adf.restart.file
            del self.adf.restart.file

        if "beckegrid" in self.adf and "atomdepquality" in self.adf.beckegrid:
            del self.adf.beckegrid.atomdepquality
            self.warnings.append('The key "beckegrid -> atomdepquality" has been removed from the input. See the new key "BeckeGrid -> QualityPerRegion".')

        if "zlmfit" in self.adf and "atomdepquality" in self.adf.zlmfit:
            del self.adf.zlmfit.atomdepquality
            self.warnings.append('The key "ZlmFit -> atomdepquality" has been removed from the input. See the new key "ZlmFit -> QualityPerRegion".')

        if "rihartreefock" in self.adf and "atomdepquality" in self.adf.rihartreefock:
            del self.adf.rihartreefock.atomdepquality
            self.warnings.append('The key "RIHartreeFock -> atomdepquality" has been removed from the input. See the new key "RIHartreeFock -> QualityPerRegion".')

    def _handle_basis(self):
        new_basis = []
        if "basis" in self.adf:
            per_atom = []
            for line in self.adf.basis._1:
                if "createoutput" in line.lower():
                    pass
                elif line.lower().startswith("type") or line.lower().startswith("core") or line.lower().startswith("fittype") or line.lower().startswith("path"):
                    new_basis.append(line.replace("/atomicdata/", "/atomicdata/ADF/"))
                elif "redfit" in line.lower():
                    self.warnings.append(f'Could not automatically convert the following line in the Basis key: {line}. "RedFit" should now be "PolTDDFT? See ADF user doc.')
                elif len(line.split()) == 2:
                    element, file = line.split()
                    per_atom.append(f"PerAtomType Symbol={element} File={file.replace('/atomicdata/', '/atomicdata/ADF/')}")
                else:
                    self.warnings.append(f'Could not automatically convert the following line in the Basis key: "{line}". See the ADF manual on the "Basis" key')

            self.adf.basis._1 = new_basis

            if per_atom:
                self.adf.basis._1 += per_atom

    def _handle_geometry_block(self):
        if "geometry" in self.adf:
            if "transitionstate" in self.adf.geometry:
                self.sett.input.ams.Task = "TransitionStateSearch"
            elif "irc" in self.adf.geometry:
                self.sett.input.ams.Task = "IRC"
            elif "lineartransit" in self.adf.geometry:
                self.sett.input.ams.Task = "PESScan"
                self.warnings.append("LinearTransit are now handled by AMS ('Task => PESScan'). See the 'PesScan' block in the AMS driver manual.")
            elif "sp" in self.adf.geometry:
                self.sett.input.ams.Task = "SinglePoint"
            elif "frequencies" in self.adf.geometry:
                self.sett.input.ams.Task = "SinglePoint"
                self.sett.input.ams.Properties.NormalModes = "Yes"
                self.sett.input.ams.NormalModes.Hessian = "Numerical"
            elif "vibron" in self.adf.geometry:
                self.warnings.append("Vibron module has been removed. For resonance Raman application see the AMS driver manual.")
            else:
                self.sett.input.ams.Task = "GeometryOptimization"

            if "converge" in self.adf.geometry and "grad" in self.adf.geometry.converge:
                self.sett.input.ams.GeometryOptimization.Convergence.Gradients = self.adf.geometry.converge.grad
            if "converge" in self.adf.geometry and "e" in self.adf.geometry.converge:
                self.sett.input.ams.GeometryOptimization.Convergence.Energy = self.adf.geometry.converge.e

            if "inithess" in self.adf.geometry:
                self.sett.input.ams.GeometryOptimization.InitialHessian.File = self.adf.geometry.inithess
                self.sett.input.ams.GeometryOptimization.InitialHessian.Type = "FromFile"

            if "mbh" in self.adf.geometry:
                self.warnings.append('Mobile Block Hessian (MBH) was not automatically converted. See "NormalModes -> Displacements Block" in the AMS driver manual.')

            del self.adf.geometry


# if __name__ == "__main__":
#     main_converter()
