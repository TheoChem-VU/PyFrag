from plams import *
import argparse as ag


__all__ = ['PyFragJob']

def ReturnArgues():
   parser = ag.ArgumentParser(description='Print user defined values')
   parser.add_argument("--ircpath", type=str, action='append', nargs='*', help='IRC coordinate file')
   parser.add_argument("--irct21", type=str, action='append', nargs='*', help='IRC coordinate file')
   parser.add_argument("--lt", type=str, action='append', nargs='*', help='LT coordinate file')
   parser.add_argument("--fragment", type=int, action='append', nargs='*', help='atom number for each fragment')
   parser.add_argument("--strain", type=float, action='append', nargs='*', help='print strain energy')
   parser.add_argument("--VDD", type=int, action='append',  nargs='*', help='print VDD charges')
   parser.add_argument("--hirshfeld", type=str, action='append',  nargs='*', help='print hirshfeld charges')
   parser.add_argument("--bondlength", type=foat, action='append', nargs='*', help='print bond length')
   parser.add_argument("--angle", type=float, action='append', nargs='*', help='print angle')
   parser.add_argument("--irrepOI", type=str, nargs='*', action='append', help='print OI energy for point group symmetry irrep')
   parser.add_argument("--population", type=str, nargs='*',action='append', help='print population for fragment orbital')
   parser.add_argument("--overlap", type=str, nargs='*',action='append', help='print overlap between two fragment orbitals')
   parser.add_argument("--orbitalenergy", type=str, nargs='*',action='append', help='print orbital energy')
   parser.add_argument("--adfinput", type=str, nargs='*',action='append', help='adfinput parameter set')
   parser.add_argument("--adfinputfile", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
   return vars(parser.parse_args())

def CollectInput(argues):
   inputKeys = {}
   for key, val in argues.items():
      if val != None:
         inputValue = []
         if key == 'overlap':
            for term in val:
               if len(term) == 4:
                  inputValue.append(({'type':term[1],'frag':term[0]},{'type':term[3],'frag':term[2]}))
               else:
                  inputValue.append(({'type':'INDEX','frag':term[1],'irrep':term[0],'index':term[2]},{'type':'INDEX','frag':term[3],'irrep':term[0],'index':term[4]}))
            inputKeys[key] = inputValue

         elif key == 'population':
            for term in val:
               if len(term) == 2:
                  inputValue.append({'type':term[1],'frag':term[0]})
               else:
                  inputValue.append({'type':'INDEX','frag':term[1],'irrep':term[0],'index':term[2]})
            inputKeys[key] = inputValue

         elif key == 'orbitalenergy':
            for term in val:
               if len(term) == 2:
                  inputValue.append({'type':term[1],'frag':term[0]})
               else:
                  inputValue.append({'type':'INDEX','frag':term[1],'irrep':term[0],'index':term[2]})
            inputKeys[key] = inputValue

         elif key == 'irrepOI':
            inputKeys[key] = [{'irrep':term[0]} for term in val]

         elif key == 'ircpath':
            inputKeys['coordFile'] = {'ircpath': val[0][0]}

         elif key == 'lt':
            inputKeys['coordFile'] = {'lt': val[0][0]}

         elif key == 'irct21':
            if len(val) == 2:
               inputKeys['coordFile'] = ({'irct21two': (val[0][0], val[1][0])})
            else:
               inputKeys['coordFile'] = ({'irct21': val[0][0]})

         elif key == 'fragment':
            inputKeys[key] = {'frag'+str(i+1):a  for i, a in enumerate(val)}

         elif key == 'strain':
            inputKeys[key] = {'frag'+str(i+1):a[0]  for i, a in enumerate(val)}

         elif key == 'VDD':
            inputKeys[key] = {'atomList': val[0]}

         elif key == 'hirshfeld':
            inputKeys[key] = [{'frag': term[0]} for term in val]

         elif key == 'bondlength':
            inputKeys[key] = [{'bondDef': [term[0], term[1]], 'oriVal': term[2]} for term in val]

         elif key == 'angle':
            inputKeys[key] = [{'angleDef': [term[0], term[1], term[2]], 'oriVal': term[3]} for term in val]

         elif key == 'adfinput':
            adfinputList   = [(term.split('=')) for term in val[0]]
            inputKeys[key] = ['sett.input.'+adfkey+'="'+keyval+'"' for adfkey, keyval in adfinputList]

         elif key == 'adfinputfile':
            f = open(val[0][0])
            adfinputLine   = [(line.split('=')) for line in f.readlines()]
            inputKeys[key] = ['sett.input.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

         else:
            inputKeys[key] = [term for term in val]

   return inputKeys


def readIRCPath(f, tag, offset):
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
   # preparations: get geometry info, atomic symbols in the right order
   f      = KFFile(fileName)
   nAtoms = f.read('Geometry', 'nr of atoms')
   aAtoms = f.read('Geometry', 'fragment and atomtype index')[nAtoms:]
   xAtoms = str(f.read('Geometry', 'atomtype')).split()
   oAtoms = f.read('Geometry', 'atom order index')[nAtoms:]
   sAtoms = [xAtoms[order-1] for numb, order in sorted(zip(oAtoms, aAtoms))]
   #sAtoms = [xAtoms[aAtoms[order-1]-1] for order in f.read('Geometry', 'atom order index')[nAtoms:]]
   bwdIRC = readIRCPath(f, 'IRC_Backward', 3*nAtoms)
   #append forward and backward coordinates
   if (len(bwdIRC) == 0) and (fileName2 != None): bwdIRC = readIRCPath(KFFile(fileName2), 'IRC_Backward', 3*nAtoms)
   fwdIRC = readIRCPath(f, 'IRC_Forward', 3*nAtoms)
   if (len(fwdIRC) == 0) and (fileName2 != None): fwdIRC = readIRCPath(KFFile(fileName2), 'IRC_Forward', 3*nAtoms)
   fwdIRC.reverse()
   #transition state geometry
   cenIRC = [f.read('IRC', 'xyz')]
   return [[[s, x, y, z] for s, x, y, z in zip(sAtoms, xyzBlock[0::3], xyzBlock[1::3], xyzBlock[2::3])] for xyzBlock in fwdIRC + cenIRC + bwdIRC]

def ReadLTt21(fileName):
   # read LT coordinates which is similar to IRC
   f      = KFFile(fileName)
   nAtoms = f.read('Geometry', 'nr of atoms')
   aAtoms = f.read('Geometry', 'fragment and atomtype index')[nAtoms:]
   xAtoms = str(f.read('Geometry', 'atomtype')).split()
   oAtoms = f.read('Geometry', 'atom order index')[nAtoms:]
   sAtoms = [xAtoms[order-1] for numb, order in sorted(zip(oAtoms, aAtoms))]
#   sAtoms = [xAtoms[aAtoms[order-1]-1] for order in f.read('Geometry', 'atom order index')[nAtoms:]]
   tmpList = f.read('LT', 'xyz')
   matrixLT = [tmpList[i:i+3*nAtoms] for i in xrange(0, tmpList, 3*nAtoms)]
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
   ircRawList = filter(None, ircRaw)
   ircFile.close()
   return ircRawList

def GetIRCFragmentList(ircStructures, fragDefinition):
   """
   #ircStructures  = from ParseIRCFile
   #fragDefinition = {"Frag1":[1,2,4], "Frag1":[3,5,6]}
   #final result will look like {'frag1':atom coordinate block, 'frag2': atom coordinate block ....}
   """
   ircList = []
   nAtoms  = sum([len(fragList) for fragList in fragDefinition.values()])
   for coordBlock in ircStructures:
      if (nAtoms != len(coordBlock)): raise RuntimeError('nAtoms in fragment definition does not match length if IRC coordinates\n')
      # loop over IRC points
      ircList.append(dict([(fragTag, Molecule()) for fragTag in fragDefinition.keys()]))
      for fragTag in fragDefinition.keys():
         # loop over fragments
         for iAtom in fragDefinition[fragTag]:
            # grab individual atoms from block according to current fragment definition
            ircList[-1][fragTag].add_atom(Atom(symbol=coordBlock[iAtom-1][0],
                                               coords=tuple([float(xyz) for xyz in coordBlock[iAtom-1][1:4]])))

   return ircList

def GetFragmentList(IRCFile, fragDefinitions):
   """
   #ircStructures  = from ParseIRCFile
   #fragDefinition = {"Frag1":[1,2,4], "Frag1":[3,5,6]}
   #final result will look like {'frag1':atom coordinate block, 'frag2': atom coordinate block ....}
   """
   ircStructures = ParseIRCFile(IRCFile)
   List = []
   # nAtoms  = sum([len(fragList) for fragList in fragDefinition.values()])
   # for coordBlock in ircStructures:
   for coordBlock,  fragDefinition in zip(ircStructures, fragDefinitions):
      ircList = []
      # if (nAtoms != len(coordBlock)): raise RuntimeError('nAtoms in fragment definition does not match length if IRC coordinates\n')
      # loop over IRC points
      ircList.append(dict([(fragTag, Molecule()) for fragTag in fragDefinition.keys()]))
      for fragTag in fragDefinition.keys():
         # loop over fragments
         for iAtom in fragDefinition[fragTag]:
            # grab individual atoms from block according to current fragment definition
            ircList[-1][fragTag].add_atom(Atom(symbol=coordBlock[iAtom-1][0],
                                               coords=tuple([float(xyz) for xyz in coordBlock[iAtom-1][1:4]])))
      complexDict = {'complex' : Molecule()}
      for i in range(len(coordBlock)):
         complexDict['complex'].add_atom(Atom(symbol=coordBlock[int(i)][0],
                                            coords=tuple([float(xyz) for xyz in coordBlock[int(i)][1:4]])))

         ircList[-1].update(complexDict)
      List.append(ircList[0])

   return List

def GetOutputTable(data):
# converge multiple values situation like {'overlap': [1.55, 1.99]} into {'overlap_1': 1.55, 'overlap_2':1.99}
   outputTable = {}
   for key, val in data.items():
      if type(val) == list and len(val) != 1:
         for i in range(len(val)):
            outputTable[key+'_'+str(i+1)] = val[i]
      elif type(val) == list and len(val) == 1:
         outputTable[key] = val[0]
      else:
         outputTable[key] = val
   return outputTable

def writeKey(file, value, pform=r'%7.5f', ljustwidth=16): # TODO: move to standalone
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
   energyfile  = open('pyfrag'+fileName+'.txt', "w")
   headerlist  = sorted(tableValues[0])
   writeKey(energyfile, headerlist) # TODO: adapt as parts of this are in standalone
   for entry in tableValues:
      sortedEntry = [entry[i] for i in headerlist]
      writeKey(energyfile, sortedEntry)
   energyfile.close()

def WriteFailFiles(failStructures, fileName):
   structureFile = open('pyfragfailed'+fileName+'.xyz', "w")
   for structure in failStructures:
      keys = structure.keys()
      structureFile.write(' '+keys[0]+' ')
      structureFile.write('\n')
      for atom in structure.values():
         for coordinate in atom:
            for term in coordinate:
               structureFile.write(' '+ term + ' ')
            structureFile.write('\n')
      structureFile.write('\n')
   structureFile.close()

def PrintTable(cellList, widthlist, bar):
#  Headwidthlist = [len(_) for _ in headersList]
#  widthlist     = [max(len(str(valuesList[_])), Headwidthlist[_]) for _ in range(len(valuesList))] # TODO: should be maximum width of the entire columns
   if bar:
      line = '-'*(sum(widthlist)+4*len(widthlist)+6)
      print ('\n', line)
   for i, entry in enumerate(cellList):
      print ('  '+str(entry).ljust(widthlist[i])+'  ') # TODO: consider number formatting here (if entry is not a string)
   print ('')
   if bar: print (line)
   # TODO: this should return a string (including "\n" if more than one line)

def PyFragDriver(inputKeys, fragmentSettings, complexSettings):
   #main pyfrag driver used for fragment and complex calculation.
   #read coordinates from IRC or LT t21 file. Other choice is xyz file generated from other tools.
   for key, val in inputKeys['coordFile'].items(): # TODO: rename coordFile to pathCoordFile
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
      outputData['StrainTotal']  = 0
      complexMolecule     = Molecule()      #each molecule is a subject of plams class Molecule()
      ircTag              = '.'+str(ircIndex+1).zfill(5)
      # for fragTag in ircFrags.keys():
      for fragTag in list(ircFrags.keys()):
         success = True
         jobFrag = ADFJob(molecule=ircFrags[fragTag], settings=fragmentSettings, name=fragTag+ircTag)
         jobFrag.run()
         if jobFrag.check():
            #provide path of fragment t21 file to final fragment analysis calculation
            exec ("complexSettings.input.Fragments." + fragTag + "=" + '"' + jobFrag.results._kfpath() + '"')
            #get strain and total strain which is the energy difference between current and previous geometry.
            outputData[fragTag + 'Strain'] = jobFrag.results.readkf('Energy', 'Bond Energy') - inputKeys['strain'][fragTag]
            outputData['StrainTotal'] += outputData[fragTag + 'Strain']
            #reorganize new complex from fragments by appending fragment label. Beware the atomic orders maybe changed
            for atom in ircFrags[fragTag]:
               atom.fragment = fragTag
               complexMolecule.add_atom(atom)
            ircFrags.pop(fragTag)
         else:
            failCases.append(ircIndex)
            success = False
            break
            # raise RuntimeError("Abnormal termination of fragment calculation, check "+fragTag+ircTag)
      if success:

         jobComplex = ADFJob(molecule=complexMolecule, settings=complexSettings, name=complexMolecule.get_formula()+ircTag)
         jobComplex.run()
         if jobComplex.check():
            successCases.append(ircIndex)
            #collect all data and put it in a list
            outputData['#IRC'] = str(ircIndex+1)
            #collect all data that need to be printed
            pyfragResult = PyFragResult(jobComplex.results, inputKeys) # TODO: See __init__ of PyFragResult below
            #convert multiple value into a dict
            outputdata   = pyfragResult.GetOutputData(complexMolecule, outputData, inputKeys)
            #collect updated informaiton of each point calculation and print it on screen
            firstIndex    =  successCases.pop(0)
            if ircIndex   == firstIndex:
               headerList = sorted(outputdata.keys())
               a = headerList.pop(headerList.index('#IRC'))
               b = headerList.pop(headerList.index('EnergyTotal'))
               headerList = [a, b] + headerList
            valuesList = [str(outputdata[i]) for i in headerList] # TODO: add formatting
            widthlist  = [max(len(str(valuesList[_])), len(str(headerList[_]))) for _ in range(len(valuesList))]
            PrintTable(headerList, widthlist, False)
            PrintTable(valuesList, widthlist, False)
            resultsTable.append(outputdata)
         else:
            success = False
      if not success:
         failStructures.append({str(ircIndex): ircStructures[ircIndex]})
         # Write faulty IRCpoints for later recovery
         if len(resultsTable) > 0:
            outputdata = {key: 'None' for key in resultsTable[0].keys() }
            outputdata.update({'#IRC': str(ircIndex+1)})
            resultsTable.append(outputdata)
   if len(resultsTable) == 0:
      raise RuntimeError("Calculations for all points failed, please check your input settings")
   return resultsTable,  str(complexMolecule.get_formula()), failStructures # return this as result only (only construct it here but use it outside)


class PyFragResult:
   def __init__(self, complexResult, inputKeys): # __init__(self, complexJob, inputKeys)
      # TODO: use complexJob here instead of complexResult and read in class members only if
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

      for key in inputKeys.keys():
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
            self.irrepOrbNum          = complexResult.readkf('Symmetry', 'norb')
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
      if len(self.irrepType) == 1:
         irreporbNum = [self.irrepOrbNum]
      else:
         irreporbNum = self.irrepOrbNum
      faIrrepone  = [[irrep for i in range(number)] for irrep, number in zip(self.irrepType, irreporbNum)]
      return  [irrep for sublist in self.faIrrepone for irrep in sublist]

   def GetOrbNum(self):
      # GetOrbNumbers including frozen core orbitals, this is necessary to read population, list like [3,4,5,12,13,14,15,23,24,26]
      #core orbital number corresponding to each irrep of complex symmetry
      coreOrbNum           = self.complexResult.readkf('Symmetry', 'ncbs')
      orbNumbers = []
      orbSum = 0
      for nrShell, nrCore in zip(self.irrepOrbNum, coreOrbNum):
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

   def GetOrbitalIndex(self, orbDescriptor):
      # orbDescriptor = {'type' = "HOMO/LUMO/INDEX", 'frag'='#frag', 'irrep'='irrepname', 'index'=i}
      fragOrbnum = self.GetFragNum(orbDescriptor['frag'])     #get fragment number
      if orbDescriptor['type'] == 'HOMO':
         index = max(xrange(len(self.orbEnergy)), key = lambda x: self.orbEnergy[x] if (self.orbFragment[x] == fragOrbnum) and self.orbOccupation[x] != 0 else -1.0E+100)
      elif orbDescriptor['type'] == 'LUMO':
         index = min(xrange(len(self.orbEnergy)), key = lambda x: self.orbEnergy[x] if (self.orbFragment[x] == fragOrbnum) and self.orbOccupation[x] == 0 else +1.0E+100)
      elif orbDescriptor['type'] == 'INDEX':
         for i in range(len(self.orbEnergy)):
            if self.orbFragment[i] == fragOrbnum  and  self.fragIrrep[i] == orbDescriptor['irrep'] and self.fragOrb[i] == int(orbDescriptor['index']):
               index = i
               break
      return index

   def ReadOverlap(self, index_1, index_2):
      #orbital numbers according to the symmetry of the complex
      faOrb = self.complexResult.readkf('SFOs', 'isfo')
      maxIndex = max(faOrb[index_1], faOrb[index_2])
      minIndex = min(faOrb[index_1], faOrb[index_2])
      index = maxIndex * (maxIndex - 1) / 2 + minIndex - 1
      if self.GetFaIrrep([index_1]) == self.GetFaIrrep([index_2]):
         self.overlap_matrix = self.complexResult.readkf(self.GetFaIrrep([index_1]), 'S-CoreSFO')
         return abs(self.overlap_matrix[index])
      else:
         return 0

   def ReadFragorbEnergy(self, index):
      return self.complexResult.readkf('Ftyp '+str(self.orbFragment[index])+self.fragIrrep[index], 'eps')[self.fragOrb[index]-1]

   def ReadIrrepOI(self, irrep):
      irrepOI              = [self.complexResult.readkf('Energy', 'Orb.Int. '+irrep)  for irrep in self.irrepType]
      fitCoefficient       = self.Int / sum(irrepOI)
      return fitCoefficient*self.complexResult.readkf('Energy', 'Orb.Int. '+irrep)

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
      outputData['EnergyTotal']  = Units.convert(self.Int, 'hartree', 'kcal/mol') + outputData['StrainTotal']
      outputData['Steric']  = outputData['Pauli'] + outputData['Elstat']
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
                     print ('bondlength', complexMolecule[atoms[0]].distance_to(complexMolecule[atoms[1]]) )
                  else:
                     value.append(Units.convert(complexMolecule[atoms[0]].distance_to(complexMolecule[atoms[1]]), 'bohr', 'angstrom') - od['oriVal'])
            outputData[key] = value

         elif key == 'angle':
            for od in val:
               atoms = self.GetAtomNum(inputKeys['fragment'], od['angleDef'])
               value.append(Units.convert((complexMolecule[atoms[0]].angle(complexMolecule[atoms[1]], complexMolecule[atoms[2]])), 'rad', 'deg') - od['oriVal'])
            outputData[key] = value
      return GetOutputTable(outputData)


# class PyFragJob:
#    def __init__(self, fragmentSettings, complexSettings, inputArgues):
#       self.inputArgues    = CollectInput(inputArgues)
#       self.fragmentSettings  = fragmentSettings
#       self.complexSettings   = complexSettings

#    def run(self):
#       tableValue, fileName, failStructures = PyFragDriver(self.inputArgues, self.fragmentSettings, self.complexSettings)
#       WriteTable(tableValue, fileName)
#       if failStructures is not None:
#          WriteFailFiles(failStructures, fileName)

class PyFragJob:
   def __init__(self, fragmentSettings, complexSettings, inputArgues, others = None):
#      if others['irct21'] == None and  others['lt'] == None and  others['ircpath'] == None:
      others_1 = dict((k,v) for k, v in others.items() if v is not None)
      others_1['irct21'] = [[inputArgues]]
#         others['irct21'] = [[inputArgues['ircpath']]]    # provide irc t21 from irc calculation
#      elif others['strain'] == None:
#         others.pop('strain')
#         others['strain'] = [[inputArgues['strain']]]    # provide strain from r1 r2 calculation
#      elif others['bondlength'][0][2] == 0:
#         others.pop('bondlength')
#         others['bondlength'] = [[inputArgues['bondlength']]]    # provide strain from r1 r2 calculation

      self.inputArgues    = CollectInput(others_1)
      self.fragmentSettings  = fragmentSettings
      self.complexSettings   = complexSettings

   def run(self):
      tableValue, fileName, failStructures = PyFragDriver(self.inputArgues, self.fragmentSettings, self.complexSettings)
      WriteTable(tableValue, fileName)
#      if failStructures is not None:
#         WriteFailFiles(failStructures, fileName)



