#from scm.plams import *
from plams import *
import argparse as ag
from  PyFragModules  import PyFragDriver, WriteTable, WriteFailFiles
"""
Pyfrag 3
Authors: Xiaobo Sun; Thomas Soini

This program has the following functionalities:
1: Reads in a series of Linear Transit or IRC structures (coordinate files or t21).
2: For each point PyFrag generates single point calculations for the individual fragments (based on user defined fragments)
   and the whole complex system.
3: In the following the corresponding ADF calculations are conducted.
4: The program will generate a text file containing the decomposition energies plus other, user defined, values such as the strain energy.

Example use:
startpython PyFrag.py  --ircpath structuresIRC_CH3N.irc --fragment 1 3 4 --fragment 2 5 --strain 0 --strain 0 --adfinput basis.type=DZ

For the earlier version (PyFrag 2.0) please see http://www.few.vu.nl/~wolters/pyfrag/
"""

parser = ag.ArgumentParser(description='Print user defined values')
parser.add_argument("--ircpath", type=str, action='append', nargs='*', help='IRC coordinate file')
parser.add_argument("--irct21", type=str, action='append', nargs='*', help='IRC coordinate file')
parser.add_argument("--lt", type=str, action='append', nargs='*', help='LT coordinate file')
parser.add_argument("--fragment", type=int, action='append', nargs='*', help='atom number for each fragment')
parser.add_argument("--strain", type=float, action='append', nargs='*', help='print strain energy')
parser.add_argument("--VDD", type=int, action='append',  nargs='*', help='print VDD charges')
parser.add_argument("--hirshfeld", type=str, action='append',  nargs='*', help='print hirshfeld charges')
parser.add_argument("--bondlength", type=float, action='append', nargs='*', help='print bond length')
parser.add_argument("--angle", type=float, action='append', nargs='*', help='print angle')
parser.add_argument("--irrepOI", type=str, nargs='*', action='append', help='print OI energy for point group symmetry irrep')
parser.add_argument("--population", type=str, nargs='*',action='append', help='print population for fragment orbital')
parser.add_argument("--overlap", type=str, nargs='*',action='append', help='print overlap between two fragment orbitals')
parser.add_argument("--orbitalenergy", type=str, nargs='*',action='append', help='print orbital energy')
# parser.add_argument("--adfinput", type=str, nargs='*',action='append', help='adfinput parameter set')
parser.add_argument("--adfinputfile", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--fragment1_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--fragment2_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--complex_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--restartjob", type=str, action='append', nargs='*', help='restart directory of file, i.e. plams.xyzuvw')
parser.add_argument("--name", type=str, action='append', nargs='*', help='provide name for the file generated by plams')
inputKeys = {'jobstate':None, 'filename':None}
for key, val in vars(parser.parse_args()).items():
   if val != None:
      inputValue = []
      if key == 'overlap':
         for term in val:
            if len(term) == 4:
               inputValue.append(({'type':term[1],'frag':term[0]},{'type':term[3],'frag':term[2]}))
            else:
               inputValue.append(({'type':'INDEX','frag':term[1],'irrep':term[0],'index':term[2]},{'type':'INDEX','frag':term[4],'irrep':term[3],'index':term[5]}))
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
         for term in val:
            if len(term) == 2:
               inputValue.append(({'bondDef': [term[0], term[1]], 'oriVal': 0}))
            else:
               inputValue.append(({'bondDef': [term[0], term[1]], 'oriVal': term[2]}))
         inputKeys[key] = inputValue

      elif key == 'angle':
         for term in val:
            if len(term) == 3:
               inputValue.append(({'angleDef': [term[0], term[1], term[2]], 'oriVal': 0}))
            else:
               inputValue.append(({'angleDef': [term[0], term[1], term[2]], 'oriVal': term[3]}))
         inputKeys[key] = inputValue

      # elif key == 'adfinput':
      #    adfinputList   = [(term.split('=')) for term in val[0]]
      #    adfComplex     = ['sett.input.'+adfkey+'="'+keyval+'"' for adfkey, keyval in adfinputList]
      #    adfFrag1       = ['sett1.input.'+adfkey+'="'+keyval+'"' for adfkey, keyval in adfinputList]
      #    adfFrag2       = ['sett2.input.'+adfkey+'="'+keyval+'"' for adfkey, keyval in adfinputList]
      #    inputKeys[key] = adfComplex + adfFrag1 + adfFrag2

      elif key == 'adfinputfile':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         adfGeneral = ['settings.input.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfFrag1 = ['settings_Frag1.input.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfFrag2 = ['settings_Frag2.input.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfComplex = ['settings_Fa.input.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         inputKeys[key] = adfGeneral + adfFrag1 + adfFrag2 + adfComplex

      elif key == 'fragment1_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_Frag1.input.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

      elif key == 'fragment2_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_Frag2.input.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

      elif key == 'complex_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_Fa.input.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]


      elif key == 'restartjob':
         inputKeys['jobstate'] = val[0][0]

      elif key == 'name':
         inputKeys['filename'] = val[0][0]

      else:
         inputKeys[key] = [term for term in val]

init(folder=inputKeys['filename'])


settings       = Settings()
settings_Frag1 = Settings()
settings_Frag2 = Settings()
settings_Fa    = Settings()


for key, val in list(inputKeys.items()):
   if key == 'adfinputfile':
      for option in inputKeys['adfinputfile']:
         exec(option)
for key, val in list(inputKeys.items()):
   if key == 'fragment1_EXTRA':
      for option in inputKeys['fragment1_EXTRA']:
         exec(option)
for key, val in list(inputKeys.items()):
   if key == 'fragment2_EXTRA':
      for option in inputKeys['fragment2_EXTRA']:
         exec(option)
for key, val in list(inputKeys.items()):
   if key == 'complex_EXTRA':
      for option in inputKeys['complex_EXTRA']:
         exec(option)

tableValue, fileName, failStructures = PyFragDriver(inputKeys, settings_Frag1, settings_Frag2, settings_Fa)

WriteTable(tableValue, fileName)
if failStructures is not None:
   WriteFailFiles(failStructures, fileName)


finish()

