import argparse as ag
import os ; os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
from qmworks import Settings, templates, run, molkit
from noodles import gather
from qmworks.packages.SCM import dftb, adf, pyfrag
from qmworks.components import Distance, select_max
from qmworks.packages.PyFragModules import GetFragmentList
#import plams
from qmworks import plams

parser = ag.ArgumentParser(description='Print user defined values')
parser.add_argument("--ircpath", type=str, action='append', nargs='*', help='IRC coordinate file')
parser.add_argument("--r1path", type=str, action='append', nargs='*', help='IRC coordinate file')
parser.add_argument("--r2path", type=str, action='append', nargs='*', help='IRC coordinate file')
parser.add_argument("--rcpath", type=str, action='append', nargs='*', help='IRC coordinate file')
parser.add_argument("--ppath", type=str, action='append', nargs='*', help='IRC coordinate file')
parser.add_argument("--tspath", type=str, action='append', nargs='*', help='IRC coordinate file')
parser.add_argument("--irct21", type=str, action='append', nargs='*', help='IRC coordinate file')
parser.add_argument("--lt", type=str, action='append', nargs='*', help='LT coordinate file')
parser.add_argument("--fragment", type=int, action='append', nargs='*', help='atom number for each fragment')
parser.add_argument("--r1fragment", type=int, action='append', nargs='*', help='atom number for each fragment')
parser.add_argument("--r2fragment", type=int, action='append', nargs='*', help='atom number for each fragment')
parser.add_argument("--rcfragment", type=int, action='append', nargs='*', help='atom number for each fragment')
parser.add_argument("--pfragment", type=int, action='append', nargs='*', help='atom number for each fragment')
parser.add_argument("--tsfragment", type=int, action='append', nargs='*', help='atom number for each fragment')
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
parser.add_argument("--R1_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--R2_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--RC_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--TS_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--P_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--IRC_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--fragment1_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--fragment2_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')
parser.add_argument("--complex_EXTRA", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')

hartree2kcal = 627.5095


inputKeys = {}
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

      elif key == 'r1path':
         inputKeys[key] = {'r1path': val[0][0]}

      elif key == 'r2path':
         inputKeys[key] = {'r2path': val[0][0]}

      elif key == 'rcpath':
         inputKeys[key] = {'rcpath': val[0][0]}

      elif key == 'ppath':
         inputKeys[key] = {'ppath': val[0][0]}

      elif key == 'tspath':
         inputKeys[key] = {'tspath': val[0][0]}

      elif key == 'lt':
         inputKeys['coordFile'] = {'lt': val[0][0]}

      elif key == 'irct21':
         if len(val) == 2:
            inputKeys['coordFile'] = ({'irct21two': (val[0][0], val[1][0])})
         else:
            inputKeys['coordFile'] = ({'irct21': val[0][0]})

      elif key == 'fragment':
         inputKeys[key] = {'frag'+str(i+1):a  for i, a in enumerate(val)}

      elif key == 'r1fragment':
         inputKeys[key] = {'frag'+str(i+1):a  for i, a in enumerate(val)}

      elif key == 'r2fragment':
         inputKeys[key] = {'frag'+str(i+1):a  for i, a in enumerate(val)}

      elif key == 'rcfragment':
         inputKeys[key] = {'frag'+str(i+1):a  for i, a in enumerate(val)}

      elif key == 'pfragment':
         inputKeys[key] = {'frag'+str(i+1):a  for i, a in enumerate(val)}

      elif key == 'tsfragment':
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
      #    inputKeys[key] = ['settings.specific.adf.'+adfkey+'="'+keyval+'"' for adfkey, keyval in adfinputList]

      elif key == 'adfinputfile':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         adfGeneral = ['settings.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfR1 = ['settings_R1.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfR2 = ['settings_R2.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfRC = ['settings_RC.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfTS = ['settings_TS.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfP = ['settings_P.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfIRC = ['settings_IRC.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfFrag1 = ['settings_Frag1.specific.fragment1.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfFrag2 = ['settings_Frag2.specific.fragment2.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         adfComplex = ['settings_Fa.specific.complex.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         inputKeys[key] = adfGeneral + adfR1 + adfR2 + adfRC + adfTS + adfIRC + adfP + adfFrag1 + adfFrag2 + adfComplex

      elif key == 'R1_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_R1.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

      elif key == 'R2_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_R2.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

      elif key == 'RC_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_RC.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

      elif key == 'TS_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_TS.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

      elif key == 'P_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_P.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

      elif key == 'IRC_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_IRC.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

      elif key == 'fragment1_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_Frag1.specific.fragment1.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

      elif key == 'fragment2_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_Frag2.specific.fragment2.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

      elif key == 'complex_EXTRA':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=',1)) for line in f.readlines()]
         inputKeys[key] = ['settings_Fa.specific.complex.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]

      else:
         inputKeys[key] = [term for term in val]

settings       = Settings()
settings_R1    = Settings()
settings_R2    = Settings()
settings_RC    = Settings()
settings_TS    = Settings()
settings_P     = Settings()
settings_IRC   = Settings()
settings_Frag1 = Settings()
settings_Frag2 = Settings()
settings_Fa    = Settings()

for key, val in list(inputKeys.items()):
   if key == 'adfinputfile':
      for option in inputKeys['adfinputfile']:
         exec(option)
for key, val in list(inputKeys.items()):
   if key == 'R1_EXTRA':
      for option in inputKeys['R1_EXTRA']:
         exec(option)
for key, val in list(inputKeys.items()):
   if key == 'R2_EXTRA':
      for option in inputKeys['R2_EXTRA']:
         exec(option)
for key, val in list(inputKeys.items()):
   if key == 'RC_EXTRA':
      for option in inputKeys['RC_EXTRA']:
         exec(option)
for key, val in list(inputKeys.items()):
   if key == 'TS_EXTRA':
      for option in inputKeys['TS_EXTRA']:
         exec(option)
for key, val in list(inputKeys.items()):
   if key == 'P_EXTRA':
      for option in inputKeys['P_EXTRA']:
         exec(option)
for key, val in list(inputKeys.items()):
   if key == 'IRC_EXTRA':
      for option in inputKeys['IRC_EXTRA']:
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

job_list = []


for key, val in list(inputKeys.items()):
   if key == 'r1path':
      ircFrags = GetFragmentList(val['r1path'], [inputKeys['r1fragment']])[0]
      r1_mol = ircFrags['frag1']
      r1 =      adf(templates.geometry.overlay(settings_R1), r1_mol, job_name="r1")
      job_list.append(gather(r1))
   elif key == 'r2path':
      ircFrags = GetFragmentList(val['r2path'], [inputKeys['r2fragment']])[0]
      r2_mol = ircFrags['frag1']
      r2 =      adf(templates.geometry.overlay(settings_R2), r2_mol, job_name="r2")
      job_list.append(gather(r2))
   elif key == 'rcpath':
      ircFrags = GetFragmentList(val['rcpath'], [inputKeys['rcfragment']])[0]
      rc_mol = ircFrags['frag1']
      rc =      adf(templates.geometry.overlay(settings_RC), rc_mol, job_name="rc")
      job_list.append(gather(rc))
   elif key == 'ppath':
      ircFrags = GetFragmentList(val['ppath'], [inputKeys['pfragment']])[0]
      p_mol = ircFrags['frag1']
      p =      adf(templates.geometry.overlay(settings_P), p_mol, job_name="p")
      job_list.append(gather(p))
   elif key == 'tspath':
      ircFrags = GetFragmentList(val['tspath'], [inputKeys['tsfragment']])[0]
      ts_mol = ircFrags['frag1']
      ts =      adf(templates.ts.overlay(settings_TS), ts_mol, job_name="ts")
      job_list.append(gather(ts))


irc = adf(templates.irc.overlay(settings_IRC), ts.molecule, job_name="irc")
pyfrag = pyfrag(templates.frag1.overlay(settings_Frag1), settings_2 = templates.frag2.overlay(settings_Frag2), settings_3 = templates.fa.overlay(settings_Fa), inputArgues = irc.kf.path , others =  vars(parser.parse_args()), job_name="pyfrag" )

job_list.append(gather(irc, pyfrag))
# Finalize and draw workflow
wf = gather(*job_list)

# Actual execution of the jobs
results = run(wf, n_processes=4)
