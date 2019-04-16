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
parser.add_argument("--adfinput", type=str, nargs='*',action='append', help='adfinput parameter set')
parser.add_argument("--adfinputfile", type=str, nargs='*',action='append', help='a file containing adfinput parameters set')


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
         inputKeys[key] = [{'bondDef': [term[0], term[1]], 'oriVal': term[2]} for term in val]

      elif key == 'angle':
         inputKeys[key] = [{'angleDef': [term[0], term[1], term[2]], 'oriVal': term[3]} for term in val]

      elif key == 'adfinput':
         adfinputList   = [(term.split('=')) for term in val[0]]
         inputKeys[key] = ['settings.specific.adf.'+adfkey+'="'+keyval+'"' for adfkey, keyval in adfinputList]

      elif key == 'adfinputfile':
         f = open(val[0][0])
         adfinputLine   = [(line.split('=')) for line in f.readlines()]
         inputKeys[key] = ['settings.specific.adf.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         inputKeys['fragmentinputfile'] = ['settings.specific.fragment.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
         inputKeys['complexinputfile'] = ['settings.specific.complex.'+adfkey+'="'+keyval.strip('\n')+'"' for adfkey, keyval in adfinputLine]
      else:
         inputKeys[key] = [term for term in val]

settings = Settings()

for key, val in list(inputKeys.items()):
   if key == 'adfinput':
      for option in inputKeys['adfinput']:
         exec(option)
   elif key == 'adfinputfile':
      for option in inputKeys['adfinputfile']:
         exec(option)
   elif key == 'fragmentinputfile':
      for option in inputKeys['fragmentinputfile']:
         exec(option)
   elif key == 'complexinputfile':
      for option in inputKeys['complexinputfile']:
         exec(option)

job_list = []


for key, val in list(inputKeys.items()):
   if key == 'r1path':
      ircFrags = GetFragmentList(val['r1path'], [inputKeys['r1fragment']])[0]
      r1_mol = ircFrags['frag1']
      r1 =      adf(templates.geometry.overlay(settings), r1_mol, job_name="r1")
      job_list.append(gather(r1))
   elif key == 'r2path':
      ircFrags = GetFragmentList(val['r2path'], [inputKeys['r2fragment']])[0]
      r2_mol = ircFrags['frag1']
      r2 =      adf(templates.geometry.overlay(settings), r2_mol, job_name="r2")
      job_list.append(gather(r2))
   elif key == 'rcpath':
      ircFrags = GetFragmentList(val['rcpath'], [inputKeys['rcfragment']])[0]
      rc_mol = ircFrags['frag1']
      rc =      adf(templates.geometry.overlay(settings), rc_mol, job_name="rc")
      job_list.append(gather(rc))
   elif key == 'ppath':
      ircFrags = GetFragmentList(val['ppath'], [inputKeys['pfragment']])[0]
      p_mol = ircFrags['frag1']
      p =      adf(templates.geometry.overlay(settings), p_mol, job_name="p")
      job_list.append(gather(p))
   elif key == 'tspath':
      ircFrags = GetFragmentList(val['tspath'], [inputKeys['tsfragment']])[0]
      ts_mol = ircFrags['frag1']
      ts =      adf(templates.ts.overlay(settings), ts_mol, job_name="ts")
      job_list.append(gather(ts))


irc = adf(templates.irc.overlay(settings), ts.molecule, job_name="irc")

#newInput={}
#newInput['ircpath'] = irc.kf.path
#newInput['strain'] = irc.kf.path
#newInput['bondlength'] = irc.kf.path

pyfrag = pyfrag(templates.frag.overlay(settings), settings_2 = templates.fa.overlay(settings), inputArgues = irc.kf.path , others =  vars(parser.parse_args()), job_name="pyfrag" )
job_list.append(gather(irc, pyfrag))
# Finalize and draw workflow
wf = gather(*job_list)

# Actual execution of the jobs
results = run(wf, n_processes=4)
