import re, os, copy
from plams import *


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
   energyfile         = open(fileName+'orbitalenergy.txt', "w")
   headerlist_all     = sorted([key for key, val in tableValues[0].items()])

   writeKey(energyfile, headerlist_all)
   for entry in tableValues:
      sortedEntry = [entry[i] for i in headerlist_all]
      writeKey(energyfile, sortedEntry)
   energyfile.close()


def GetreactionFiles(fpath, fragment):
   '''This function is used to recursively traverse a directory and find
   all .t21 files contained within. It is used as a generator.'''
   for dir, subdirs, files in os.walk(fpath):
      subdirs.sort(reverse = False)
      for subdir in subdirs:
         for p in GetreactionFiles(subdir,fragment):
            yield p
      files.sort(reverse = False)
      for file in files:
         if file.startswith(fragment) and file.endswith(".t21"):
            yield os.path.join(dir, file)


def OrbitalEnergy(inputKeys):
   energylist    = []
   for key, val in inputKeys['fragment'].items():
      file_generator   = GetreactionFiles(inputKeys['coordFile']['irct21'],val[0])
      for i, file in enumerate(file_generator):
         kf            = KFFile(file)
         orbitalEnergy = PyFragResult(kf)
         energy = orbitalEnergy.GetOutputData({}, inputKeys)
         fileExtension = os.path.splitext(file)[0]
         fileNum       = os.path.splitext(fileExtension)[-1]
         energy['#IRC']  = fileNum
         energylist.append(energy)
   return energylist

def ConvertList(obj):
   #single number in adf t21 is number fommat which is need to convert list
   if type(obj) == list:
      return obj
   else:
      return [obj]


class PyFragResult:
   def __init__(self, complexResult):

      #number of orbitals for each symmetry for complex
      self.irrepOrbNumber       = ConvertList(complexResult.read('Symmetry', 'norb'))

      #total number of orbitals

      self.totalOrb             =  sum(self.irrepOrbNumber)

      #irrep label for symmetry of complex
      self.irrepType            = str(complexResult.read('Symmetry', 'symlab')).split()

      #energy for each orbital
      self.orbEnergy            = sum([ConvertList(complexResult.read(self.irrepType[i], 'eps_A'))  for i in range(len(self.irrepType))], []) + \
                                  sum([ConvertList(complexResult.read(self.irrepType[i], 'eps_B'))  for i in range(len(self.irrepType))], [])
      #occupation of each orbitals which is either 0 or 1
      self.orbOccupation        = sum([ConvertList(complexResult.read(self.irrepType[i], 'froc_A'))  for i in range(len(self.irrepType))], []) + \
                                  sum([ConvertList(complexResult.read(self.irrepType[i], 'froc_B'))  for i in range(len(self.irrepType))], [])

      # expand the orb nubber according to irrepType(0,1,2,3,4,    0,1,2,3,4)
      self.irrepOrblist         = [num+1  for i in range(len(self.irrepOrbNumber)) for num in range(self.irrepOrbNumber[i])] + \
                                  [num+1  for i in range(len(self.irrepOrbNumber)) for num in range(self.irrepOrbNumber[i])]

      # [a1, a1, b1, b1, b2, b2]
      self.irrepTypelist        = sum( [[irrep for i in range(number)] for irrep, number in zip(self.irrepType, self.irrepOrbNumber)], []) + \
                                  sum( [[irrep for i in range(number)] for irrep, number in zip(self.irrepType, self.irrepOrbNumber)], [])


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
      orbIndex = 0
      if self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'HOMO':
         orbIndex = sorted(range(len(self.orbEnergy)), key = lambda x: self.orbEnergy[x] if  self.orbOccupation[x] != 0 else -1.0E+100, reverse=True)[-int(self.GetFrontIndex(orbDescriptor['type'])['num'])]
      elif self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'LUMO':
         orbIndex = sorted(range(len(self.orbEnergy)), key = lambda x: self.orbEnergy[x] if  self.orbOccupation[x] == 0 else +1.0E+100)[int(self.GetFrontIndex(orbDescriptor['type'])['num'])]
      elif self.GetFrontIndex(orbDescriptor['type'])['holu'] == 'INDEX':
         for i in range(len(self.orbEnergy)):
            if  self.irrepTypelist[i] == orbDescriptor['irrep'] and self.irrepOrblist[i] == int(orbDescriptor['index']):
               orbIndex = i
               break
      return orbIndex


   def GetOutputData(self, outputData, inputKeys):
      #collect user defined data
      for key, val in list(inputKeys.items()):
         value = []

         if key == 'orbitalenergy':
            for od in val:
               if od['type'] != 'INDEX':
                  outputData[od['type']]                             = Units.convert(self.orbEnergy[self.GetOrbitalIndex(od)], 'hartree', 'eV')
               else:
                  outputData[od['irrep'] + '_' + od['index'] + '_A'] = Units.convert(self.orbEnergy[self.GetOrbitalIndex(od)], 'hartree', 'eV')
                  outputData[od['irrep'] + '_' + od['index'] + '_B'] = Units.convert(self.orbEnergy[self.GetOrbitalIndex(od) + self.totalOrb], 'hartree', 'eV')
      return outputData
