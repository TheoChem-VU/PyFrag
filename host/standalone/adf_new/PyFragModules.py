from scm.plams import AMSJob, KFFile, Molecule, Atom, Units, load_all, Settings, AMSResults, log
# from scm.libbase import InputParser  # type: ignore
# from scm.plams.core.functions import parse_heredoc
import re
import os
from typing import List, Dict, Optional, Tuple, Any
import logging
# from plams import *


def settings_from_inputfile(inputfile: str) -> Settings:
    """
    Returns a Settings object from an inputfile

    Adapted from the AMSJob.from_inputfile method to read in the settings from a file and return a Settings object
    Reason for adapting was because the AMSJob.from_inputfile method does remove the ams block from the settings object
    """
    pre_calc_job: AMSJob = AMSJob.from_inputfile(inputfile)

    # This does not include the "System" block; it does include the "Task" and "Engine" block
    settings = pre_calc_job.settings

    # # First make sure that the inputfile is read in correctly (parsing the heredoc)
    # with open(inputfile, 'r') as f:
    #     inp_file = parse_heredoc(f.read(), 'eor')

    # # Then read in the "System" block from the inputfile and add it to the settings object
    # unnested_settings = InputParser().to_settings("ams", inp_file)
    # settings.input["ams"].update(unnested_settings["ams"])

    # Work around for the fact that the ams block is not included in the settings object
    # The above code apparently does not allow for overriding the molecule
    if pre_calc_job.molecule is not None:
        molecule: Molecule = pre_calc_job.molecule[""]  # type: ignore
        charge = molecule.properties.charge
        settings.input.ams.System.Charge = charge

    return settings


def HandleRestart(foldername) -> Optional[str]:
    """
    Copied mostly from the plams routine on restarting jobs
    https://github.com/SCM-NV/PLAMS/blob/master/scripts/plams
    Input:
       foldername (str): the foldername to restart from
    """
    restart_backup = None
    if os.path.isdir(foldername):
        foldername = foldername.rstrip('/')
        if os.listdir(foldername):
            restart_backup = f"{foldername}.res"
            n = 1
            while os.path.exists(restart_backup):
                n += 1
                restart_backup = f"{foldername}.res{str(n)}"
            os.rename(foldername, restart_backup)
            print(f'RESTART: Moving {foldername} to {restart_backup} and restarting from it')
    else:
        restart_backup = None
        print('RESTART: The folder specified for restart does not exist.')
    return restart_backup


def ReadIRCPath(f, tag, offset):
    # split all data into different block each of which represent one molecular structure in IRC
    if f.read(tag, 'PathStatus').strip() == 'DONE' or f.read(tag, 'PathStatus').strip() == 'EXEC':
        nEntry = f.read(tag, 'CurrentPoint') * offset
        tmpList = f.read(tag, 'xyz')
        return [tmpList[i:i+offset] for i in range(0, nEntry, offset)]
    return []


def GetAtom(kf) -> Tuple[int, List[str]]:
    # preparations: get geometry info, atomic symbols in the right order
    nAtoms = kf.read('Geometry', 'nr of atoms')
    aAtoms = kf.read('Geometry', 'fragment and atomtype index')[nAtoms:]
    xAtoms = str(kf.read('Geometry', 'atomtype')).split()
    oAtoms = kf.read('Geometry', 'atom order index')[nAtoms:]
    sAtoms = [xAtoms[order-1] for numb, order in sorted(zip(oAtoms, aAtoms))]
    # sAtoms = [xAtoms[aAtoms[order-1]-1] for order in f.read('Geometry', 'atom order index')[nAtoms:]]
    return nAtoms, sAtoms


def ReadIRCt21(fileName, fileName2=None):
    f = KFFile(fileName)
    nAtoms, sAtoms = GetAtom(f)
    bwdIRC = ReadIRCPath(f, 'IRC_Backward', 3*nAtoms)
    # append forward and backward coordinates
    if (len(bwdIRC) == 0) and (fileName2 is not None):
        bwdIRC = ReadIRCPath(KFFile(fileName2), 'IRC_Backward', 3*nAtoms)
    fwdIRC = ReadIRCPath(f, 'IRC_Forward', 3*nAtoms)
    if (len(fwdIRC) == 0) and (fileName2 is not None):
        fwdIRC = ReadIRCPath(KFFile(fileName2), 'IRC_Forward', 3*nAtoms)
    fwdIRC.reverse()
    # transition state geometry
    cenIRC = [f.read('IRC', 'xyz')]
    return [[[s, x, y, z] for s, x, y, z in zip(sAtoms, xyzBlock[0::3], xyzBlock[1::3], xyzBlock[2::3])] for xyzBlock in fwdIRC + cenIRC + bwdIRC]


def ReadLTt21(fileName):
    # read LT coordinates which is similar to IRC
    f = KFFile(fileName)
    nAtoms, sAtoms = GetAtom(f)
    tmpList = f.read('LT', 'xyz')
    matrixLT = [tmpList[i:i+3*nAtoms] for i in range(0, len(tmpList), 3*nAtoms)]
    return [[[s, x, y, z] for s, x, y, z in zip(sAtoms, xyzBlock[0::3], xyzBlock[1::3], xyzBlock[2::3])] for xyzBlock in matrixLT]


def ParseIRCFile(ircCoordFile):
    # read xyz file like amv file exported from adfinput
    ircFile = open(str(ircCoordFile))
    ircRaw = [[]]
    for line in ircFile:
        llist = line.split()
        lllen = len(llist)
        if lllen == 4:
            # append coordinate
            ircRaw[-1].append(llist)
        else:
            # initiate new coordinate block
            ircRaw.append([])
    ircRawList = [_f for _f in ircRaw if _f]
    ircFile.close()
    return ircRawList


def GetIRCFragmentList(ircStructures, fragDefinition):
    """
    #ircStructures  = from ParseIRCFile
    #fragDefinition = {"Frag1":[1,2,4], "Frag1":[3,5,6]}
    #final result will look like {'frag1':atom coordinate block, 'frag2': atom coordinate block ....}
    """
    ircList = []
    nAtoms = sum([len(fragList) for fragList in list(fragDefinition.values())])
    for coordBlock in ircStructures:
        if (nAtoms != len(coordBlock)):
            raise RuntimeError('nAtoms in fragment definition does not match length if IRC coordinates\n')
        # loop over IRC points
        ircList.append(dict([(fragTag, Molecule()) for fragTag in list(fragDefinition.keys())]))
        for fragTag in list(fragDefinition.keys()):
            # loop over fragments
            for iAtom in fragDefinition[fragTag]:
                # grab individual atoms from block according to current fragment definition
                ircList[-1][fragTag].add_atom(Atom(symbol=coordBlock[iAtom-1][0],
                                                   coords=tuple([float(xyz) for xyz in coordBlock[iAtom-1][1:4]])))
    return ircList


def GetOutputTable(data):
    # converge multiple values situation like {'overlap': [1.55, 1.99]} into {'overlap_1': 1.55, 'overlap_2':1.99}
    outputTable = {}
    for key, val in list(data.items()):
        if isinstance(val, list) and len(val) != 1:
            for i in range(len(val)):
                outputTable[key+'_'+str(i+1)] = val[i]
        elif isinstance(val, list) and len(val) == 1:
            outputTable[key] = val[0]
        else:
            outputTable[key] = val
    return outputTable


def writeKey(file, value, pform=r'%7.5f', ljustwidth=16):
    # write all data into a file.Keep 7 digits and 5 decimals and the width of each entry is 16
    for val in value:
        if val is None:
            file.write(str.ljust('  ---', ljustwidth))
        else:
            if isinstance(val, float):
                file.write(str.ljust(pform % (val), ljustwidth))
            else:
                file.write(str.ljust(str(val), ljustwidth))
    file.write('\n')


def WriteTable(tableValues, fileName):
    energyfile = open(f"pyfrag_{fileName}.txt", "w")
    headerlist_all = sorted(tableValues[0])
    # check if bondlength exist
    bondlist = [e for e in headerlist_all if e in ("bondlength")]
    if len(bondlist) == 1:
        # order the print list
        headerlist_select = [e for e in headerlist_all if e not in ("#IRC", "bondlength", "EnergyTotal", "Int", "Elstat", "Pauli", "OI", "Disp", "StrainTotal", "frag1Strain", "frag2Strain")]
        headerlist = ["#IRC", "bondlength", "EnergyTotal", "Int", "Elstat", "Pauli", "OI", "Disp", "StrainTotal", "frag1Strain", "frag2Strain"] + headerlist_select
    else:
        headerlist_select = [e for e in headerlist_all if e not in ("#IRC", "bondlength_1", "EnergyTotal", "Int", "Elstat", "Pauli", "OI", "Disp", "StrainTotal", "frag1Strain", "frag2Strain")]
        headerlist = ["#IRC", "bondlength_1", "EnergyTotal", "Int", "Elstat", "Pauli", "OI", "Disp", "StrainTotal", "frag1Strain", "frag2Strain"] + headerlist_select
    writeKey(energyfile, headerlist)
    for entry in tableValues:
        sortedEntry = [entry[i] for i in headerlist]
        writeKey(energyfile, sortedEntry)
    energyfile.close()


def WriteFailFiles(failStructures, fileName):
    structureFile = open(f"pyfragfailed_{fileName}.xyz", "w")
    for structure in failStructures:
        keys = list(structure.keys())
        structureFile.write(' '+keys[0]+' ')
        structureFile.write('\n')
        for atom in list(structure.values()):
            for coordinate in atom:
                for term in coordinate:
                    structureFile.write(' ' + term + ' ')
                structureFile.write('\n')
        structureFile.write('\n')
    structureFile.close()


def PrintTable(cellList, widthlist, bar):
    if bar:
        line = '-'*(sum(widthlist)+4*len(widthlist)+6)
        print('\n', line)
    for i, entry in enumerate(cellList):
        print('  '+str(entry).ljust(widthlist[i])+'  ', end=" ")
    print('')
    if bar:
        print(line)


def CleanUpCalculationFolder(job: AMSJob):
    """
    Removes unnecessary files from the calculation folder
    See https://www.scm.com/doc/plams/components/results.html#scm.plams.core.results.Results._clean for more info
    """
    r: AMSResults = job.results
    if r.ok():
        mol = r.get_main_molecule()
        r._clean(["-", "$JN.err", "$JN.run", "CreateAtoms.out", "t12.rel", "output.xyz"])
        for atom in set(mol.atoms):
            r._clean(["-", f"t21.*.{atom.symbol}"])
        job.pickle()  # this will update the .dill file which is used to restart the job and extract results when using plams


def PyFragDriver(inputKeys, frag1Settings, frag2Settings, complexSettings):
    # main pyfrag driver used for fragment and complex calculation.
    # read coordinates from IRC or LT t21 file. Other choice is xyz file generated from other tools.
    if inputKeys['jobstate'] is not None:
        load_all(inputKeys['jobstate'])

    for key, val in list(inputKeys['coordFile'].items()):
        if key == 'irct21':
            ircStructures = ReadIRCt21(val)
            exec('complexSettings.input.UNITS.length="Bohr"')
        elif key == 'irct21two':
            ircStructures = ReadIRCt21(val[0], val[1])
            exec('complexSettings.input.UNITS.length="Bohr"')
        elif key == 'lt':
            ircStructures = ReadLTt21(val)
            exec('complexSettings.input.UNITS.length="Bohr"')
        else:
            ircStructures = ParseIRCFile(val)

    resultsTable = []
    failCases = []
    successCases = []
    failStructures = []

    for ircIndex, ircFrags in enumerate(GetIRCFragmentList(ircStructures, inputKeys['fragment'])):
        logging.log(level=logging.INFO, msg=f"Starting calculations for IRC point {ircIndex+1}/{len(ircStructures)}")
        logging.log(level=logging.DEBUG, msg="\n".join(str(ircFrags[fragtag]) for fragtag in sorted(list(ircFrags.keys()))))
        outputData = {}
        outputData['StrainTotal'] = 0
        ircTag = '.'+str(ircIndex+1).zfill(5)
        for fragTag in sorted(list(ircFrags.keys())):
            success = True
            if fragTag == 'frag1':
                fragmentSettings = frag1Settings
            else:
                fragmentSettings = frag2Settings
            for coorKey, coorVal in list(inputKeys['coordFile'].items()):
                if coorKey != 'ircpath':
                    exec('fragmentSettings.input.UNITS.length="Bohr"')
            if fragTag == 'frag1':
                jobFrag1 = AMSJob(molecule=ircFrags[fragTag], settings=fragmentSettings, name=fragTag+ircTag)
                jobFrag1.run()
                frag1Molecule = ircFrags[fragTag]

                # check if the calculation is successful and log error message if not
                if not jobFrag1.results.ok():
                    logging.log(level=logging.CRITICAL, msg=f"Fragment calculation for {fragTag} failed, please check your input settings")

                outputData[fragTag + 'Strain'] = jobFrag1.results.get_energy(unit='kcal/mol') - inputKeys['strain'][fragTag]
                outputData['StrainTotal'] += outputData[fragTag + 'Strain']
                CleanUpCalculationFolder(jobFrag1)
            else:
                jobFrag2 = AMSJob(molecule=ircFrags[fragTag], settings=fragmentSettings, name=fragTag+ircTag)
                jobFrag2.run()
                frag2Molecule = ircFrags[fragTag]

                # check if the calculation is successful and log error message if not
                if not jobFrag2.results.ok():
                    logging.log(level=logging.CRITICAL, msg=f"Fragment calculation for {fragTag} failed, please check your input settings")

                outputData[fragTag + 'Strain'] = jobFrag2.results.get_energy(unit='kcal/mol') - inputKeys['strain'][fragTag]
                outputData['StrainTotal'] += outputData[fragTag + 'Strain']
                CleanUpCalculationFolder(jobFrag2)

            # disable the result check because ADF print a lot of error message
            # if jobFrag.check():
            if True:
                ircFrags.pop(fragTag)
            else:
                failCases.append(ircIndex)
                success = False
                break
        if success:  # succes if always true, needs to be fixed

            for at in frag1Molecule:
                at.properties.suffix = 'adf.f=frag1'

            for at in frag2Molecule:
                at.properties.suffix = 'adf.f=frag2'

            complexMolecule = frag1Molecule+frag2Molecule
            complexSettings.input.adf.fragments.frag1 = (jobFrag1, 'adf')
            complexSettings.input.adf.fragments.frag2 = (jobFrag2, 'adf')
            jobComplex = AMSJob(molecule=complexMolecule, settings=complexSettings, name='complex'+ircTag)
            logging.log(level=logging.INFO, msg=f"Running complex {ircIndex+1}")
            jobComplex.run()

            if not jobComplex.results.ok():
                logging.log(level=logging.CRITICAL, msg="Complex calculationfailed, please check your input settings")

            # disable the result check because ADF print a lot of useless message
            if True:
                # if jobComplex.check():
                successCases.append(ircIndex)
                # collect all data and put it in a list
                outputData['#IRC'] = str(ircIndex+1)
                # collect all data that need to be printed
                pyfragResult = PyFragResult(jobComplex.results, inputKeys)
                # convert multiple value into a dict
                outputLine = pyfragResult.GetOutputData(complexMolecule, outputData, inputKeys)
                # collect updated informaiton of each point calculation and print it on screen
                firstIndex = successCases.pop(0)
                if ircIndex == firstIndex:
                    headerList = sorted(outputLine.keys())
                    a = headerList.pop(headerList.index('#IRC'))
                    b = headerList.pop(headerList.index('EnergyTotal'))
                    headerList = [a, b] + headerList
                valuesList = [str(outputLine[i]) for i in headerList]
                widthlist = [max(len(str(valuesList[_])), len(str(headerList[_]))) for _ in range(len(valuesList))]
                PrintTable(headerList, widthlist, False)
                PrintTable(valuesList, widthlist, False)
                resultsTable.append(outputLine)
                CleanUpCalculationFolder(jobComplex)
                logging.log(level=logging.INFO, msg=f"IRC point {ircIndex+1}/{len(ircStructures)} finished")
            else:
                success = False
        if not success:
            failStructures.append({str(ircIndex): ircStructures[ircIndex]})
            # Write faulty IRCpoints for later recovery
            if len(resultsTable) > 0:
                outputLine = {key: 'None' for key in list(resultsTable[0].keys())}
                outputLine.update({'#IRC': str(ircIndex+1)})
                resultsTable.append(outputLine)
    if len(resultsTable) == 0:
        raise RuntimeError("Calculations for all points failed, please check your input settings")
    return resultsTable, inputKeys['filename'], failStructures  # return this as result only (only construct it here but use it outside)


class PyFragResult:
    def __init__(self, complexResult, inputKeys):  # __init__(self, complexJob, inputKeys)
        # 1. needed for output requested by user 2. complexJob.check passes
        self.complexResult = complexResult
        # Pauli energy
        self.Pauli = complexResult.readrkf('Energy', 'Pauli Total', file='adf')
        # Electrostatic energy
        self.Elstat = complexResult.readrkf('Energy', 'Electrostatic Interaction', file='adf')
        # total OI which is usually not equal to the sum of all irrep OI
        self.OI = complexResult.readrkf('Energy', 'Orb.Int. Total', file='adf')
        # energy of total complex which is the sum of Pauli, Elstat and OI
        self.Int = complexResult.readrkf('Energy', 'Bond Energy', file='adf')
        # Dispersion Energy
        self.Disp = complexResult.readrkf('Energy', 'Dispersion Energy', file='adf')
        # irrep label for symmetry of complex
        self.irrepType = str(complexResult.readrkf('Symmetry', 'symlab', file='adf')).split()

        for key in list(inputKeys.keys()):
            if key == 'overlap' or key == 'population' or key == 'orbitalenergy' or key == 'irrepOI':
                # orbital numbers according to the symmetry of each fragment and the orbitals belonging to the same symmetry in different fragments
                self.fragOrb = complexResult.readrkf('SFOs', 'ifo', file='adf')
                # symmetry for each orbital of fragments
                self.fragIrrep = str(complexResult.readrkf('SFOs', 'subspecies', file='adf')).split()
                # the fragment label for each orbital
                self.orbFragment = complexResult.readrkf('SFOs', 'fragment', file='adf')
                # energy for each orbital
                self.orbEnergy = complexResult.readrkf('SFOs', 'energy', file='adf')
                # occupation of each orbitals which is either 0 or 2
                self.orbOccupation = complexResult.readrkf('SFOs', 'occupation', file='adf')
                # number of orbitals for each symmetry for complex
                self.irrepOrbNumber = complexResult.readrkf('Symmetry', 'norb', file='adf')
                # number of core orbitals for each symmetry for complex
                self.coreOrbNumber = complexResult.readrkf('Symmetry', 'ncbs', file='adf')

    def ConvertList(self, obj) -> List[Any]:
        # single number in adf t21 is number fommat which is need to convert list
        if isinstance(obj, list):
            return obj
        else:
            return [obj]

    def GetFaIrrep(self):
        # append complex irrep label to each orbital, if symmetry is A, convert self.irrepOrbNum which is float type into list
        irreporbNum = self.ConvertList(self.irrepOrbNumber)
        faIrrepone = [[irrep for i in range(number)] for irrep, number in zip(self.irrepType, irreporbNum)]
        return [irrep for sublist in faIrrepone for irrep in sublist]

    def GetOrbNum(self):
        # GetOrbNumbers including frozen core orbitals, this is necessary to read population, list like [3,4,5,12,13,14,15,23,24,26]
        # core orbital number corresponding to each irrep of complex symmetry
        coreOrbNum = self.ConvertList(self.coreOrbNumber)
        irrepOrbNum = self.ConvertList(self.irrepOrbNumber)
        orbNumbers = []
        orbSum = 0
        for nrShell, nrCore in zip(irrepOrbNum, coreOrbNum):
            orbSum += (nrShell + nrCore)
            orbNumbers.extend(list(range(orbSum - nrShell + 1, orbSum + 1)))
        return orbNumbers

    def GetFragOrbNum(self):
        # GetOrbNumbers including frozen core orbitals, this is necessary to read population, list like [3,4,5,12,13,14,15,23,24,26]
        # core orbital number corresponding to each irrep of complex symmetry
        coreOrbNum = self.ConvertList(self.coreOrbNumber)
        irrepOrbNum = self.ConvertList(self.irrepOrbNumber)
        orbNumbers = []
        for nrShell, nrCore in zip(irrepOrbNum, coreOrbNum):
            orbNumbers.extend(range(nrCore + 1, nrShell + nrCore + 1))
        return orbNumbers

    def GetAtomNum(self, fragmentList, atoms) -> List[int]:
        # change atom number in old presentation of a molecule into atom number in new presentation that formed by assembling fragments
        atomList = [atomNum for key in sorted(list(fragmentList.keys())) for atomNum in list(fragmentList[key])]
        return [atomList.index(i)+1 for i in atoms]

    def GetFragNum(self, frag: str) -> int:
        # change frag type like 'frag1' into number like "1" recorded in t21
        # append fragmenttype(like 1 or 2) to each orbital
        fragType = str(self.complexResult.readrkf('Geometry', 'fragmenttype', file='adf')).split()
        return fragType.index(frag) + 1

    def GetFrontIndex(self, orbSign) -> Dict[str, str]:
        # convert HOMO/LUMO/HOMO-1/LUMO+1/INDEX into dict {'holu': 'HOMO', 'num': -1}
        for matchString in [r'HOMO(.*)', r'LUMO(.*)', r'INDEX']:
            matchObj = re.match(matchString, orbSign)
            if matchObj:
                holu = re.sub(r'(.[0-9]+)', "", matchObj.group())
                num = re.sub(r'([a-zA-Z]+)', "", matchObj.group())
                if num:
                    return {'holu': holu, 'num': num}
                else:
                    return {'holu': holu, 'num': 0}

    def GetOrbitalIndex(self, orbDescriptor):
        # orbDescriptor = {'type' = "HOMO/LUMO/INDEX", 'frag'='#frag', 'irrep'='irrepname', 'index'=i}
        fragOrbnum = self.GetFragNum(orbDescriptor['frag'])  # get fragment number
        orbIndex = 0
        if self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'HOMO':
            orbIndex = sorted(range(len(self.orbEnergy)), key=lambda x: self.orbEnergy[x] if (self.orbFragment[x] == fragOrbnum)
                              and self.orbOccupation[x] != 0 else -1.0E+100, reverse=True)[-int(self.GetFrontIndex(orbDescriptor['type'])['num'])]
        elif self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'LUMO':
            orbIndex = sorted(range(len(self.orbEnergy)), key=lambda x: self.orbEnergy[x] if (self.orbFragment[x] == fragOrbnum)
                              and self.orbOccupation[x] == 0 else +1.0E+100)[int(self.GetFrontIndex(orbDescriptor['type'])['num'])]
        elif self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'INDEX':
            for i in range(len(self.orbEnergy)):
                if self.orbFragment[i] == fragOrbnum and self.fragIrrep[i] == orbDescriptor['irrep'] and self.fragOrb[i] == int(orbDescriptor['index']):
                    orbIndex = i
                    break
        return orbIndex

    def ReadOverlap(self, index_1, index_2) -> float:
        # orbital numbers according to the symmetry of the complex
        faOrb = self.GetFragOrbNum()
        faIrrep = self.GetFaIrrep()
        maxIndex = max(faOrb[index_1], faOrb[index_2])
        minIndex = min(faOrb[index_1], faOrb[index_2])
        index = maxIndex * (maxIndex - 1) / 2 + minIndex - 1
        if faIrrep[index_1] == faIrrep[index_2]:
            self.overlap_matrix = self.complexResult.readrkf(faIrrep[index_1], 'S-CoreSFO', file='adf')
            return abs(self.overlap_matrix[int(index)])
        else:
            return 0

    def ReadFragorbEnergy(self, index):
        return self.complexResult.readrkf('Ftyp '+str(self.orbFragment[index])+self.fragIrrep[index], 'eps', file='adf')[self.fragOrb[index]-1]

    def CheckIrrepOI(self) -> List[Dict[str, str]]:
        """ Checks whether OI irreps are present that can be included. Returns a list of irreps that can be included."""
        irreps = [{"irrep": irrep} for irrep in self.irrepType]
        logging.log(level=logging.INFO, msg=f"Found irreps [{[irrep.values() for irrep in irreps]}] in complex that will be included in OI")
        return irreps

    def ReadIrrepOI(self, irrep) -> float:
        irrepOI = [self.complexResult.readrkf('Energy', 'Orb.Int. '+irreps, file='adf') for irreps in self.irrepType]
        fitCoefficient = self.OI / sum(irrepOI)
        return fitCoefficient*Units.convert(self.complexResult.readrkf('Energy', 'Orb.Int. '+irrep, file='adf'), 'hartree', 'kcal/mol')

    def ReadPopulation(self, index) -> float:
        orbNumbers = self.GetOrbNum()
        # populations of all orbitals
        sfoPopul = self.complexResult.readrkf('SFO popul', 'sfo_grosspop', file='adf')
        return sfoPopul[orbNumbers[index] - 1]

    def ReadVDD(self, atomList) -> List[float]:
        vddList = []
        for atom in atomList:
            vddScf = self.complexResult.readrkf('Properties', 'AtomCharge_SCF Voronoi', file='adf')[int(atom)-1]
            vddInit = self.complexResult.readrkf('Properties', 'AtomCharge_initial Voronoi', file='adf')[int(atom)-1]
            vddList.append(vddScf - vddInit)
        return vddList

    def ReadHirshfeld(self, fragment) -> float:
        valueHirshfeld = self.complexResult.readrkf('Properties', 'FragmentCharge Hirshfeld', file='adf')
        return valueHirshfeld[self.GetFragNum(fragment) - 1]

    def GetOutputData(self, complexMolecule, outputData, inputKeys):
        # collect default energy parts for activation strain analysis
        outputData['Pauli'] = Units.convert(self.Pauli, 'hartree', 'kcal/mol')
        outputData['Elstat'] = Units.convert(self.Elstat, 'hartree', 'kcal/mol')
        outputData['OI'] = Units.convert(self.OI, 'hartree', 'kcal/mol')
        outputData['Int'] = Units.convert(self.Int, 'hartree', 'kcal/mol')
        outputData['Disp'] = Units.convert(self.Disp, 'hartree', 'kcal/mol')
        outputData['EnergyTotal'] = Units.convert(self.Int, 'hartree', 'kcal/mol') + outputData['StrainTotal']

        # check for unspecified options such as irrep printing if not specified by user
        if 'irrepOI' not in inputKeys and len(self.irrepType) != 1:
            log(level=3, message="NOTE:No irreps specified for OI calculation by the user, but detecting irreps in complex.")
            inputKeys["irrepOI"] = self.CheckIrrepOI()

        # collect user defined data
        for key, val in list(inputKeys.items()):
            value = []
            if key == 'overlap':
                outputData[key] = [self.ReadOverlap(self.GetOrbitalIndex(od1), self.GetOrbitalIndex(od2)) for od1, od2 in val]

            elif key == 'population':
                outputData[key] = [self.ReadPopulation(self.GetOrbitalIndex(od)) for od in val]

            elif key == 'orbitalenergy':
                outputData[key] = [self.ReadFragorbEnergy(self.GetOrbitalIndex(od)) for od in val]

            elif key == 'irrepOI':
                outputData[key] = [self.ReadIrrepOI(od['irrep']) for od in val]

            elif key == 'hirshfeld':
                outputData[key] = [self.ReadHirshfeld(od['frag']) for od in val]

            elif key == 'VDD':
                outputData[key] = self.ReadVDD(val['atomList'])

            elif key == 'bondlength':
                for od in val:
                    atoms = self.GetAtomNum(inputKeys['fragment'], od['bondDef'])
                    # coordinate read directly from t21 is in bohr while from .amv is in amstrom
                    for coorKey, coorVal in list(inputKeys['coordFile'].items()):
                        if coorKey == 'ircpath':
                            value.append(complexMolecule[atoms[0]].distance_to(complexMolecule[atoms[1]]) - od['oriVal'])
                            # print ('bondlength', complexMolecule[atoms[0]].distance_to(complexMolecule[atoms[1]]) )
                        else:
                            value.append(Units.convert(complexMolecule[atoms[0]].distance_to(complexMolecule[atoms[1]]), 'bohr', 'angstrom') - od['oriVal'])
                outputData[key] = value

            elif key == 'angle':
                for od in val:
                    atoms = self.GetAtomNum(inputKeys['fragment'], od['angleDef'])
                    value.append(Units.convert((complexMolecule[atoms[0]].angle(complexMolecule[atoms[1]], complexMolecule[atoms[2]])), 'rad', 'deg') - od['oriVal'])
                outputData[key] = value
        return GetOutputTable(outputData)
