#from scm.plams import *
import re, copy
from plams import *
def ReadIRCPath(f, tag, offset):
   # split all data into different block each of which represent one molecular structure in IRC
   if f.read(tag, 'PathStatus').strip() == 'DONE' or f.read(tag, 'PathStatus').strip() == 'EXEC':
      nEntry  = f.read(tag, 'CurrentPoint') * offset
      tmpList = f.read(tag, 'xyz')
      return [tmpList[i:i+offset] for i in range(0, nEntry, offset)]
   return []

def GetAtom(kf):
   # preparations: get geometry info, atomic symbols in the right order
   nAtoms = kf.read('Geometry', 'nr of atoms')
   aAtoms = kf.read('Geometry', 'fragment and atomtype index')[nAtoms:]
   xAtoms = str(kf.read('Geometry', 'atomtype')).split()
   oAtoms = kf.read('Geometry', 'atom order index')[nAtoms:]
   sAtoms = [xAtoms[order-1] for numb, order in sorted(zip(oAtoms, aAtoms))]
   #sAtoms = [xAtoms[aAtoms[order-1]-1] for order in f.read('Geometry', 'atom order index')[nAtoms:]]
   return nAtoms, sAtoms

def ReadIRCt21(fileName, fileName2=None):
   f              = KFFile(fileName)
   nAtoms, sAtoms = GetAtom(f)
   bwdIRC         = ReadIRCPath(f, 'IRC_Backward', 3*nAtoms)
   #append forward and backward coordinates
   if (len(bwdIRC) == 0) and (fileName2 != None): bwdIRC = ReadIRCPath(KFFile(fileName2), 'IRC_Backward', 3*nAtoms)
   fwdIRC = ReadIRCPath(f, 'IRC_Forward', 3*nAtoms)
   if (len(fwdIRC) == 0) and (fileName2 != None): fwdIRC = ReadIRCPath(KFFile(fileName2), 'IRC_Forward', 3*nAtoms)
   fwdIRC.reverse()
   #transition state geometry
   cenIRC = [f.read('IRC', 'xyz')]
   return [[[s, x, y, z] for s, x, y, z in zip(sAtoms, xyzBlock[0::3], xyzBlock[1::3], xyzBlock[2::3])] for xyzBlock in fwdIRC + cenIRC + bwdIRC]

def ReadLTt21(fileName):
   # read LT coordinates which is similar to IRC
   f              = KFFile(fileName)
   nAtoms, sAtoms = GetAtom(f)
   tmpList        = f.read('LT', 'xyz')
   matrixLT       = [tmpList[i:i+3*nAtoms] for i in range(0, tmpList, 3*nAtoms)]
   return [[[s, x, y, z] for s, x, y, z in zip(sAtoms, xyzBlock[0::3], xyzBlock[1::3], xyzBlock[2::3])] for xyzBlock in tmpList]


def ParseIRCFile(ircCoordFile):
   # read xyz file like amv file exported from adfinput
   ircFile = open(str(ircCoordFile))
   ircRaw  = [[]]
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
   nAtoms  = sum([len(fragList) for fragList in list(fragDefinition.values())])
   for coordBlock in ircStructures:
      if (nAtoms != len(coordBlock)): raise RuntimeError('nAtoms in fragment definition does not match length if IRC coordinates\n')
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
      if type(val) == list and len(val) != 1:
         for i in range(len(val)):
            outputTable[key+'_'+str(i+1)] = val[i]
      elif type(val) == list and len(val) == 1:
         outputTable[key] = val[0]
      else:
         outputTable[key] = val
   return outputTable

def writeKey(file, value, pform=r'%7.3f', ljustwidth=16):
   # write all data into a file.Keep 7 digits and 5 decimals and the width of each entry is 16
   for val in value:
      if val is None:
         file.write(str.ljust('  ---', ljustwidth))
      else:
         if type(val) == float:
            file.write(str.ljust(pform % (val), ljustwidth))
         else:
            file.write(str.ljust(str(val), ljustwidth))
   file.write('\n')

def WriteTable(tableValues, fileName):
   energyfile         = open('pyfrag1'+fileName+'.txt', "w")
   headerlist_all     = sorted(tableValues[0])
   #order the print list
   headerlist_select  = [e for e in headerlist_all if e not in ("#IRC","bondlength","EnergyTotal","StrainTotal","Int","Elstat","Pauli","OI","Disp","int","straintotal","frag1Strain","frag2Strain")]
   headerlist         = ["#IRC","bondlength","EnergyTotal","StrainTotal","Int","Elstat","Pauli","OI","Disp"] + headerlist_select
   writeKey(energyfile, headerlist)
   for entry in tableValues:
      sortedEntry = [entry[i] for i in headerlist]
      writeKey(energyfile, sortedEntry)
   energyfile.close()

def WriteDefaultTable(tableValues, fileName):
   energyfile  = open('pyfrag2'+fileName+'.txt', "w")
   headerlist  = ["#IRC","bondlength","EnergyTotal","int","straintotal","frag1Strain","frag2Strain"]
   writeKey(energyfile, headerlist)
   for entry in tableValues:
      sortedEntry = [entry[i] for i in headerlist]
      writeKey(energyfile, sortedEntry)
   energyfile.close()


def WriteFailFiles(failStructures, fileName):
   structureFile = open('pyfragfailed'+fileName+'.xyz', "w")
   for structure in failStructures:
      keys = list(structure.keys())
      structureFile.write(' '+keys[0]+' ')
      structureFile.write('\n')
      for atom in list(structure.values()):
         for coordinate in atom:
            for term in coordinate:
               structureFile.write(' '+ term + ' ')
            structureFile.write('\n')
      structureFile.write('\n')
   structureFile.close()

def PrintTable(cellList, widthlist, bar):
   if bar:
      line = '-'*(sum(widthlist)+4*len(widthlist)+6)
      print ('\n', line)
   for i, entry in enumerate(cellList):
      print ('  '+str(entry).ljust(widthlist[i])+'  ')
   print ('')
   if bar: print (line)

def PyFragDriver(inputKeys, frag1Settings, frag2Settings, complexSettings, frag1Settings_open, frag2Settings_open, complexSettings_open):
   #main pyfrag driver used for fragment and complex calculation.
   #read coordinates from IRC or LT t21 file. Other choice is xyz file generated from other tools.
#   print (complexSettings, '\n', frag1Settings, '\n', frag2Settings)
   if inputKeys['jobstate'] is not None:
      load_all(inputKeys['jobstate'])

   for key, val in list(inputKeys['coordFile'].items()):
      if key == 'irct21':
         ircStructures = ReadIRCt21(val)
         exec ('complexSettings.input.UNITS.length="Bohr"')
      elif key == 'irct21two':
         ircStructures = ReadIRCt21(val[0],val[1])
         exec ('complexSettings.input.UNITS.length="Bohr"')
      elif key == 'lt':
         ircStructures = ReadLTt21(val)
         exec ('complexSettings.input.UNITS.length="Bohr"')
      else:
         ircStructures = ParseIRCFile(val)


   resultsTable   = []
   failCases      = []
   successCases   = []
   failStructures = []

   for ircIndex, ircFrags in enumerate(GetIRCFragmentList(ircStructures, inputKeys['fragment'])):
      outputData = {}
      outputData_open = {}
      outputData['StrainTotal']  = 0
      outputData_open['StrainTotal'] = 0
      outputData_open['straintotal'] = 0
      ircFrags_open = copy.deepcopy(ircFrags)
      complexMolecule          = Molecule()      #each molecule is a subject of plams class Molecule()
      complexMolecule_open     = Molecule()
      ircTag              = '.'+str(ircIndex+1).zfill(5)
      # for fragTag in ircFrags.keys():
      for fragTag in sorted(list(ircFrags.keys())):
         if fragTag == 'frag1':
            fragmentSettings = frag1Settings
         else:
            fragmentSettings = frag2Settings
         for coorKey, coorVal in list(inputKeys['coordFile'].items()):
            if coorKey != 'ircpath':
               exec ('fragmentSettings.input.UNITS.length="Bohr"')

         jobFrag = ADFJob(molecule=ircFrags[fragTag], settings=fragmentSettings, name=fragTag+ircTag)
         jobFrag.run()

         #provide path of fragment t21 file to final fragment analysis calculation
         if inputKeys['jobstate'] is not None:
            exec ("complexSettings.input.Fragments." + fragTag + "=" + '"' +  inputKeys['jobstate'] + "/" + fragTag + ircTag + "/" + fragTag + ircTag + ".t21" + '"')
         else:
            exec ("complexSettings.input.Fragments." + fragTag + "=" + '"' + jobFrag.results._kfpath() + '"')
         #get strain and total strain which is the energy difference between current and previous geometry.
         # outputData[fragTag + 'Strain'] = Units.convert(jobFrag.results.readkf('Energy', 'Bond Energy'), 'hartree', 'kcal/mol') - inputKeys['strain'][fragTag]
         # outputData['StrainTotal'] += outputData[fragTag + 'Strain']
         #reorganize new complex from fragments by appending fragment label. Beware the atomic orders maybe changed
         for atom in ircFrags[fragTag]:
            atom.fragment = fragTag
            complexMolecule.add_atom(atom)
         ircFrags.pop(fragTag)


      jobComplex = ADFJob(molecule=complexMolecule, settings=complexSettings, name=complexMolecule.get_formula()+ircTag)
      jobComplex.run()
      if jobComplex.check():
         successCases.append(ircIndex)
         #collect all data and put it in a list
         outputData['#IRC'] = str(ircIndex+1)
         #collect all data that need to be printed
         pyfragResult = PyFragResult(jobComplex.results, inputKeys)
         #convert multiple value into a dict
         outputLine   = pyfragResult.GetOutputData(complexMolecule, outputData, inputKeys)
         #collect updated informaiton of each point calculation and print it on screen
         firstIndex    =  successCases.pop(0)
         if ircIndex   == firstIndex:
            headerList = sorted(outputLine.keys())
            a = headerList.pop(headerList.index('#IRC'))
            b = headerList.pop(headerList.index('EnergyTotal'))
            headerList = [a, b] + headerList
         valuesList = [str(outputLine[i]) for i in headerList]
         widthlist  = [max(len(str(valuesList[_])), len(str(headerList[_]))) for _ in range(len(valuesList))]
         PrintTable(headerList, widthlist, False)
         PrintTable(valuesList, widthlist, False)


      for fragTag_open in sorted(list(ircFrags_open.keys())):
         if fragTag_open == 'frag1':
            fragmentSettings = frag1Settings_open
         else:
            fragmentSettings = frag2Settings_open
         for coorKey, coorVal in list(inputKeys['coordFile'].items()):
            if coorKey != 'ircpath':
               exec ('fragmentSettings.input.UNITS.length="Bohr"')

         jobFrag = ADFJob(molecule=ircFrags_open[fragTag_open], settings=fragmentSettings, name=fragTag_open+"open"+ircTag)
         jobFrag.run()

         #get strain and total strain which is the energy difference between current and previous geometry.
         outputData_open[fragTag_open + 'Strain'] = Units.convert(jobFrag.results.readkf('Energy', 'Bond Energy'), 'hartree', 'kcal/mol') - inputKeys['strain'][fragTag_open]
         outputData_open['straintotal'] += outputData_open[fragTag_open + 'Strain']

         #Reorganize new complex from fragments by appending fragment label. Beware the atomic orders maybe changed
         for atom in ircFrags_open[fragTag_open]:
            complexMolecule_open.add_atom(atom)
         ircFrags_open.pop(fragTag_open)

      jobComplex = ADFJob(molecule=complexMolecule_open, settings=complexSettings_open, name=complexMolecule_open.get_formula()+ "open"+ircTag)
      jobComplex.run()
      if jobComplex.check():
         #collect all data that need to be printed
         # outputData_open['EnergyTotal'] = Units.convert(jobComplex.results.readkf('Energy', 'Bond Energy'), 'hartree', 'kcal/mol') - inputKeys['straintotal']['straintotal']
         outputData_open['EnergyTotal'] = Units.convert(jobComplex.results.readkf('Energy', 'Bond Energy'), 'hartree', 'kcal/mol') - inputKeys['strain']['frag1'] - inputKeys['strain']['frag2']
         outputData_open['int'] = outputData_open['EnergyTotal'] - outputData_open['straintotal']
         outputData_open['StrainTotal'] = outputData_open['EnergyTotal'] - outputLine['Int']
         #convert multiple value into a dict
         outputLine_open   = GetOutputTable(outputData_open)
         outputLine.update(outputLine_open)
         #collect updated informaiton of each point calculation and print it on screen
         resultsTable.append(outputLine)

   if len(resultsTable) == 0:
      raise RuntimeError("Calculations for all points failed, please check your input settings")

   return resultsTable,  str(complexMolecule.get_formula()), failStructures # return this as result only (only construct it here but use it outside)


class PyFragResult:
   def __init__(self, complexResult, inputKeys): # __init__(self, complexJob, inputKeys)
      # 1. needed for output requested by user 2. complexJob.check passes
      self.complexResult        = complexResult
      #Pauli energy
      self.Pauli                = complexResult.readkf('Energy', 'Pauli Total')
      #Electrostatic energy
      self.Elstat               = complexResult.readkf('Energy', 'Electrostatic Interaction')
      #total OI which is usually not equal to the sum of all irrep OI
      self.OI                   = complexResult.readkf('Energy', 'Orb.Int. Total')
      #energy of total complex which is the sum of Pauli, Elstat and OI
      self.Int                  = complexResult.readkf('Energy', 'Bond Energy')
      #Dispersion Energy
      self.Disp                  = complexResult.readkf('Energy', 'Dispersion Energy')

      for key in list(inputKeys.keys()):
         if key == 'overlap' or key == 'population' or key == 'orbitalenergy' or key == 'irrepOI':
         #orbital numbers according to the symmetry of each fragment and the orbitals belonging to the same symmetry in different fragments
            self.fragOrb              = complexResult.readkf('SFOs', 'ifo')
            #symmetry for each orbital of fragments
            self.fragIrrep            = str(complexResult.readkf('SFOs', 'subspecies')).split()
            #the fragment label for each orbital
            self.orbFragment          = complexResult.readkf('SFOs', 'fragment')
            #energy for each orbital
            self.orbEnergy            = complexResult.readkf('SFOs', 'energy')
            #occupation of each orbitals which is either 0 or 2
            self.orbOccupation        = complexResult.readkf('SFOs', 'occupation')
            #number of orbitals for each symmetry for complex
            self.irrepOrbNumber          = complexResult.readkf('Symmetry', 'norb')
            #irrep label for symmetry of complex
            self.irrepType            = str(complexResult.readkf('Symmetry', 'symlab')).split()
            self.coreOrbNumber        = complexResult.readkf('Symmetry', 'ncbs')

   def ConvertList(self, obj):
      #single number in adf t21 is number fommat which is need to convert list
      if type(obj) == list:
         return obj
      else:
         return [obj]

   def GetFaIrrep(self):
      #append complex irrep label to each orbital, if symmetry is A, convert self.irrepOrbNum which is float type into list
      irreporbNum = self.ConvertList(self.irrepOrbNumber)
      faIrrepone  = [[irrep for i in range(number)] for irrep, number in zip(self.irrepType, irreporbNum)]
      return  [irrep for sublist in faIrrepone for irrep in sublist]

   def GetOrbNum(self):
      # GetOrbNumbers including frozen core orbitals, this is necessary to read population, list like [3,4,5,12,13,14,15,23,24,26]
      #core orbital number corresponding to each irrep of complex symmetry
      coreOrbNum           = self.ConvertList(self.coreOrbNumber)
      irrepOrbNum          = self.ConvertList(self.irrepOrbNumber)
      orbNumbers = []
      orbSum = 0
      for nrShell, nrCore in zip(irrepOrbNum, coreOrbNum):
         orbSum += (nrShell + nrCore)
         orbNumbers.extend(list(range(orbSum - nrShell + 1, orbSum + 1)))
      return orbNumbers

   def GetFragOrbNum(self):
      # GetOrbNumbers including frozen core orbitals, this is necessary to read population, list like [3,4,5,12,13,14,15,23,24,26]
      #core orbital number corresponding to each irrep of complex symmetry
      coreOrbNum           = self.ConvertList(self.coreOrbNumber)
      irrepOrbNum          = self.ConvertList(self.irrepOrbNumber)
      orbNumbers = []
      for nrShell, nrCore in zip(irrepOrbNum, coreOrbNum):
         orbNumbers.extend(range(nrCore + 1, nrShell + nrCore + 1))
      return orbNumbers

   def GetAtomNum(self, fragmentList, atoms):
      #change atom number in old presentation of a molecule into atom number in new presentation that formed by assembling fragments
      atomList = [atomNum for key in sorted(list(fragmentList.keys())) for atomNum in list(fragmentList[key])]
      return [atomList.index(i)+1 for i in atoms]

   def GetFragNum(self, frag):
      #change frag type like 'frag1' into number like "1" recorded in t21
      #append fragmenttype(like 1 or 2) to each orbital
      fragType = str(self.complexResult.readkf('Geometry', 'fragmenttype')).split()
      return fragType.index(frag) + 1

   def GetFrontIndex(self, orbSign):
     #convert HOMO/LUMO/HOMO-1/LUMO+1/INDEX into dict {'holu': 'HOMO', 'num': -1}
      for matchString in [r'HOMO(.*)', r'LUMO(.*)', r'INDEX']:
         matchObj = re.match(matchString, orbSign)
         if matchObj:
            holu = re.sub(r'(.[0-9]+)',"", matchObj.group())
            num = re.sub(r'([a-zA-Z]+)',"", matchObj.group())
            if num:
               return {'holu': holu, 'num': num}
            else:
               return {'holu': holu, 'num': 0}

   def GetOrbitalIndex(self, orbDescriptor):
      # orbDescriptor = {'type' = "HOMO/LUMO/INDEX", 'frag'='#frag', 'irrep'='irrepname', 'index'=i}
      fragOrbnum = self.GetFragNum(orbDescriptor['frag'])     #get fragment number
      orbIndex = 0
      if self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'HOMO':
         orbIndex = sorted(range(len(self.orbEnergy)), key = lambda x: self.orbEnergy[x] if (self.orbFragment[x] == fragOrbnum) and self.orbOccupation[x] != 0 else -1.0E+100, reverse=True)[-int(self.GetFrontIndex(orbDescriptor['type'])['num'])]
      elif self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'LUMO':
         orbIndex = sorted(range(len(self.orbEnergy)), key = lambda x: self.orbEnergy[x] if (self.orbFragment[x] == fragOrbnum) and self.orbOccupation[x] == 0 else +1.0E+100)[int(self.GetFrontIndex(orbDescriptor['type'])['num'])]
      elif self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'INDEX':
         for i in range(len(self.orbEnergy)):
            if self.orbFragment[i] == fragOrbnum  and  self.fragIrrep[i] == orbDescriptor['irrep'] and self.fragOrb[i] == int(orbDescriptor['index']):
               orbIndex = i
               break
      print ("orbIndex",orbIndex)
      return orbIndex

   def ReadOverlap(self, index_1, index_2):
      #orbital numbers according to the symmetry of the complex
      faOrb    = self.GetFragOrbNum()
      faIrrep  = self.GetFaIrrep()
      maxIndex = max(faOrb[index_1], faOrb[index_2])
      minIndex = min(faOrb[index_1], faOrb[index_2])
      index = maxIndex * (maxIndex - 1) / 2 + minIndex - 1
      if faIrrep[index_1] == faIrrep[index_2]:
         self.overlap_matrix = self.complexResult.readkf(faIrrep[index_1], 'S-CoreSFO')
         return abs(self.overlap_matrix[int(index)])
      else:
         return 0

   def ReadFragorbEnergy(self, index):
      return self.complexResult.readkf('Ftyp '+str(self.orbFragment[index])+self.fragIrrep[index], 'eps')[self.fragOrb[index]-1]

   def ReadIrrepOI(self, irrep):
      irrepOI              = [self.complexResult.readkf('Energy', 'Orb.Int. '+irreps)  for irreps in self.irrepType]
      fitCoefficient       = self.OI / sum(irrepOI)
      return fitCoefficient*Units.convert(self.complexResult.readkf('Energy', 'Orb.Int. '+irrep), 'hartree', 'kcal/mol')

   def ReadPopulation(self, index):
      orbNumbers =  self.GetOrbNum()
      #populations of all orbitals
      sfoPopul = self.complexResult.readkf('SFO popul', 'sfo_grosspop')
      return sfoPopul[orbNumbers[index] - 1]

   def ReadVDD(self, atomList):
      vddList = []
      for atom in atomList:
         vddScf  = self.complexResult.readkf('Properties', 'AtomCharge_SCF Voronoi')[int(atom)-1]
         vddInit = self.complexResult.readkf('Properties', 'AtomCharge_initial Voronoi')[int(atom)-1]
         vddList.append(vddScf - vddInit)
      return vddList

   def ReadHirshfeld(self, fragment):
      valueHirshfeld = self.complexResult.readkf('Properties', 'FragmentCharge Hirshfeld')
      return valueHirshfeld[self.GetFragNum(fragment) - 1]

   def GetOutputData(self, complexMolecule, outputData, inputKeys):
      #collect default energy parts for activation strain analysis
      outputData['Pauli']   = Units.convert(self.Pauli, 'hartree', 'kcal/mol')
      outputData['Elstat']  = Units.convert(self.Elstat, 'hartree', 'kcal/mol')
      outputData['OI']      = Units.convert(self.OI, 'hartree', 'kcal/mol')
      outputData['Int']     = Units.convert(self.Int, 'hartree', 'kcal/mol')
      outputData['Disp']     = Units.convert(self.Disp, 'hartree', 'kcal/mol')
      outputData['EnergyTotal']  = Units.convert(self.Int, 'hartree', 'kcal/mol') + outputData['StrainTotal']
      #collect user defined data
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
               #coordinate read directly from t21 is in bohr while from .amv is in amstrom
               for coorKey, coorVal in list(inputKeys['coordFile'].items()):
                  if coorKey == 'ircpath':
                     value.append(complexMolecule[atoms[0]].distance_to(complexMolecule[atoms[1]]) - od['oriVal'])
                     #print ('bondlength', complexMolecule[atoms[0]].distance_to(complexMolecule[atoms[1]]) )
                  else:
                     value.append(Units.convert(complexMolecule[atoms[0]].distance_to(complexMolecule[atoms[1]]), 'bohr', 'angstrom') - od['oriVal'])
            outputData[key] = value

         elif key == 'angle':
            for od in val:
               atoms = self.GetAtomNum(inputKeys['fragment'], od['angleDef'])
               value.append(Units.convert((complexMolecule[atoms[0]].angle(complexMolecule[atoms[1]], complexMolecule[atoms[2]])), 'rad', 'deg') - od['oriVal'])
            outputData[key] = value
      return GetOutputTable(outputData)

