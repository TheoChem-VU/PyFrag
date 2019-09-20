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

def GetIRCFragmentList(ircStructures):
   """
   #ircStructures  = from ParseIRCFile
   #fragDefinition = {"Frag1":[1,2,4], "Frag2":[3,5,6]}
   #final result will look like {'frag1':atom coordinate block, 'frag2': atom coordinate block ....}
   """
   ircList = []
   for coordBlock in ircStructures:
      # loop over IRC points
      ircList.append(Molecule())
      for iAtom in coordBlock:
         # grab individual atoms from block according to current fragment definition
         ircList[-1].add_atom(Atom(symbol=iAtom[0], coords=tuple([float(xyz) for xyz in iAtom[1:4]])))

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
   energyfile  = open('pyfrag'+fileName+'.txt', "w")
   headerlist  = sorted(tableValues[0])
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

def PyFragDriver(inputKeys, complexSettings):
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

   for ircIndex, ircFrags in enumerate(GetIRCFragmentList(ircStructures)):
      outputData = {}
      outputData_open = {}
      complexMolecule          = Molecule()      #each molecule is a subject of plams class Molecule()
      ircTag              = '.'+str(ircIndex+1).zfill(5)
      jobComplex = ADFJob(molecule=ircFrags, settings=complexSettings, name=ircFrags.get_formula()+ircTag)
      jobComplex.run()
      if jobComplex.check():
         successCases.append(ircIndex)
         #collect all data and put it in a list
         outputData['#IRC'] = str(ircIndex+1)
         #collect all data that need to be printed
         pyfragResult = PyFragResult(jobComplex.results, inputKeys)
         #convert multiple value into a dict
         outputLine   = pyfragResult.GetOutputData(ircFrags, outputData, inputKeys)

         #collect updated informaiton of each point calculation and print it on screen
         firstIndex    =  successCases.pop(0)
         if ircIndex   == firstIndex:
            headerList = sorted(outputLine.keys())
            a = headerList.pop(headerList.index('#IRC'))
            headerList = [a] + headerList
         valuesList = [str(outputLine[i]) for i in headerList]
         widthlist  = [max(len(str(valuesList[_])), len(str(headerList[_]))) for _ in range(len(valuesList))]
         PrintTable(headerList, widthlist, False)
         PrintTable(valuesList, widthlist, False)

         #collect all data that need to be printed
         outputData_open['EnergyTotal'] = Units.convert(jobComplex.results.readkf('Energy', 'Bond Energy'), 'hartree', 'kcal/mol')
         #convert multiple value into a dict
         outputLine_open   = GetOutputTable(outputData_open)
         outputLine.update(outputLine_open)
         #collect updated informaiton of each point calculation and print it on screen
         resultsTable.append(outputLine)

   if len(resultsTable) == 0:
      raise RuntimeError("Calculations for all points failed, please check your input settings")

   return resultsTable,  str(ircFrags.get_formula()), failStructures # return this as result only (only construct it here but use it outside)

class PyFragResult:
   def __init__(self, complexResult, inputKeys): # __init__(self, complexJob, inputKeys)
      # 1. needed for output requested by user 2. complexJob.check passes
      self.complexResult        = complexResult

   def ConvertList(self, obj):
      #single number in adf t21 is number fommat which is need to convert list
      if type(obj) == list:
         return obj
      else:
         return [obj]

   def ReadVDD(self, atomList):
      vddList = []
      for atom in atomList:
         vddScf  = self.complexResult.readkf('Properties', 'AtomCharge_SCF Voronoi')[int(atom)-1]
         vddInit = self.complexResult.readkf('Properties', 'AtomCharge_initial Voronoi')[int(atom)-1]
         vddList.append(vddScf - vddInit)
      return vddList

   def GetOutputData(self, complexMolecule, outputData, inputKeys):
      #collect default energy parts for activation strain analysis
      #collect user defined data
      for key, val in list(inputKeys.items()):
         value = []

         if key == 'VDD':
            outputData[key] = self.ReadVDD(val['atomList'])

         elif key == 'bondlength':
            for od in val:
               atoms = od['bondDef']
               #coordinate read directly from t21 is in bohr while from .amv is in amstrom
               for coorKey, coorVal in list(inputKeys['coordFile'].items()):
                  if coorKey == 'ircpath':
                     value.append(complexMolecule[int(atoms[0])].distance_to(complexMolecule[int(atoms[1])]) - od['oriVal'])
                     #print ('bondlength', complexMolecule[atoms[0]].distance_to(complexMolecule[atoms[1]]) )
                  else:
                     value.append(Units.convert(complexMolecule[int(atoms[0])].distance_to(complexMolecule[int(atoms[1])]), 'bohr', 'angstrom') - od['oriVal'])
            outputData[key] = value

         elif key == 'angle':
            for od in val:
               atoms = od['angleDef']
               value.append(Units.convert((complexMolecule[int(atoms[0])].angle(complexMolecule[int(atoms[1])], complexMolecule[int(atoms[2])])), 'rad', 'deg') - od['oriVal'])
            outputData[key] = value
      return GetOutputTable(outputData)
