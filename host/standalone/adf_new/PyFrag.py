from scm.plams import *
import os
import sys
import shutil
# from plams import *
import argparse as ag
import importlib.util
from  PyFragModules  import PyFragDriver, WriteTable, WriteFailFiles, HandleRestart

"""
Pyfrag 3
Authors: Xiaobo Sun; Thomas Soini; Siebe Lekanne Deprez

This program has the following functionalities:
1: Reads in a series of Linear Transit or IRC structures (coordinate files or t21).
2: For each point PyFrag generates single point calculations for the individual fragments (based on user defined fragments)
   and the whole complex system.
3: In the following the corresponding ADF calculations are conducted.
4: The program will generate a text file containing the decomposition energies plus other, user defined, values such as the strain energy.

Example use:
amspython PyFrag.py  --ircpath structuresIRC_CH3N.irc --fragment 1 3 4 --fragment 2 5 --strain 0 --strain 0 --adfinput basis.type=DZ

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
parser.add_argument("--adfinputfile", type=str, nargs='*',action='append', help='a file containing adfinput parameters set of ADF after 2019')
parser.add_argument("--old_adfinputfile", type=str, nargs='*',action='append', help='a file containing adfinput parameters set of ADF prior to 2019')
parser.add_argument("--fragment1_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--fragment2_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--complex_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--restartjob", type=str, action='append', nargs='*', help='restart directory of file, i.e. plams.xyzuvw')
parser.add_argument("--name", type=str, action='append', nargs='*', help='provide name for the file generated by plams')
inputKeys = {'jobstate':None, 'filename':None}
for key, val in vars(parser.parse_args()).items():
   if val != None:
      # print(f"Key = {key} \nValue = {val}") # for debugging
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

      elif (key == 'adfinputfile' or key == 'old_adfinputfile' or key == 'fragment1_EXTRA' or key == 'fragment2_EXTRA' or key == 'complex_EXTRA'):
         inputKeys[key] = val[0][0]   

      elif key == 'restartjob':
         # val[0][0] is the restart directory which may be changed with a .0x extension by the HandleRestart function
         # function is imported from the PyFragModules.py module
         inputKeys['jobstate'] = HandleRestart(val[0][0]) 

      elif key == 'name':
         inputKeys['filename'] = val[0][0]

      else:
         inputKeys[key] = [term for term in val]


init(folder=inputKeys['filename'])
workdir_path = config.default_jobmanager.workdir


# This (ugly) block of code is necessary to check if the user has provided an AMS input file or an old ADF input file
old_ADF_input = False
if "adfinputfile" in inputKeys.keys():
   settings_general = AMSJob.from_inputfile(inputKeys['adfinputfile']).settings
elif "old_adfinputfile" in inputKeys.keys():
   print("Detecting an old ADF inputfile (prior to 2019).\nAttempting to convert the input to settings WHICH MAY NOT BE SUCCESFUL...")
   sys.path.append("/".join(__file__.rsplit("/")[:-3]))  # path to the adf_to_ams_input_converter.py file
   from adf_to_ams_input_converter import main_converter
   settings_general = main_converter((inputKeys['old_adfinputfile']))
   old_ADF_input = True

settings_Frag1   = settings_general.copy() # copy is necessary to avoid linking errors when changing settings
settings_Frag2   = settings_general.copy()
settings_Complex = settings_general.copy()

for extra_input, extra_input_path in inputKeys.items():
   if extra_input == 'fragment1_EXTRA':
      if old_ADF_input:
         settings_Frag1 += main_converter(extra_input_path)
      else:
         settings_Frag1 += AMSJob.from_inputfile(extra_input_path).settings
   if extra_input == 'fragment2_EXTRA':
      if old_ADF_input:
         settings_Frag2 += main_converter(extra_input_path)
      else:
         settings_Frag2 += AMSJob.from_inputfile(extra_input_path).settings
   if extra_input == 'complex_EXTRA':
      if old_ADF_input:
         settings_Complex += main_converter(extra_input_path)
      else:
         settings_Complex += AMSJob.from_inputfile(extra_input_path).settings

tableValue, fileName, failStructures = PyFragDriver(inputKeys, settings_Frag1, settings_Frag2, settings_Complex)

WriteTable(tableValue, fileName)
if failStructures is not None:
   WriteFailFiles(failStructures, fileName)

# Remove extra files
[os.remove(os.path.join(workdir_path.rsplit('/', 1)[0], filename)) 
   for filename in ["sub", "adfinputfile", "old_adfinputfile", "complex_EXTRA", "fragment1_EXTRA", "fragment2_EXTRA"]
   if os.path.exists(os.path.join(workdir_path.rsplit('/', 1)[0], filename)) 
] 

# Remove the restart directory if it exists. Note: this may be quite a dangerous operation
if inputKeys['jobstate'] is not None:
    shutil.rmtree(inputKeys['jobstate'])
finish()

