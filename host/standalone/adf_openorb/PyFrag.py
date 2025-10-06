import argparse as ag
from  PyFragModules  import WriteTable, OrbitalEnergy

parser = ag.ArgumentParser(description='Print user defined values')
parser.add_argument("--ircpath", type=str, action='append', nargs='*', help='IRC coordinate file')
parser.add_argument("--orbitalenergy", type=str, nargs='*',action='append', help='print orbital energy')
parser.add_argument("--fragment", type=str, action='append', nargs='*', help='which fragment file')

'''
inputKeys = {'coordFile': {'irct21': '/Users/xiaobo/Desktop/pyfragopen/plams.0001'}, 'orbitalenergy': [{'type': 'HOMO'}, {'type': 'LUMO'}, {'type': 'INDEX', 'irrep': 'AA', 'index': '5'}], 'fragment': {'frag1': ['frag1open']}}
'''
inputKeys = {}
for key, val in vars(parser.parse_args()).items():
   if val != None:
      inputValue = []

      if key == 'orbitalenergy':
         for term in val:
            if len(term) == 1:
               inputValue.append({'type':term[0]})
            else:
               inputValue.append({'type':'INDEX','irrep':term[0],'index':term[1]})
         inputKeys[key] = inputValue

      elif key == 'ircpath':
         inputKeys['coordFile'] = ({'irct21': val[0][0]})

      elif key == 'fragment':
         inputKeys[key] = {'frag'+str(i+1):a  for i, a in enumerate(val)}

      else:
         inputKeys[key] = [term for term in val]


'''
[{'HOMO': -5.867204517099804, 'LUMO': -3.4971059493027004, 'AA_5_A': -270.203926788353, 'AA_5_B': -269.84360190285327, '#IRC': '.00001'}, {'HOMO': -5.867204517099804, 'LUMO': -3.4971059493027004, 'AA_5_A': -270.203926788353, 'AA_5_B': -269.84360190285327, '#IRC': '.00010'}]
'''

tableValues = OrbitalEnergy(inputKeys)
WriteTable(tableValues, inputKeys['fragment']['frag1'][0])
