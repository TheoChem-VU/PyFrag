#!/usr/bin/env python

def read_energy_t21(t21, unit='kcal/mol'):

    input = open(t21, 'r')
    gout = input.read()
    input.close()
    energyline=re.compile(r'total energy      ='+r'.*?'+ r':  kinetic energy', re.IGNORECASE|re.DOTALL).findall(gout)
    energy = energyline[0].split()
    return float(energy[3])*627.51

def bondlength_xyz(xyz, atom1, atom2):

    try:
        bond = [int(atom1), int(atom2)]
    except:
        sys.stderr.write("values for printing bondlength are not integers!")
        return None

    bondlength = math.sqrt((xyz[bond[0]-1][0]-xyz[bond[1]-1][0])**2 +\
                       (xyz[bond[0]-1][1]-xyz[bond[1]-1][1])**2 +\
                       (xyz[bond[0]-1][2]-xyz[bond[1]-1][2])**2)
    return bondlength

def shelljoin(*args):
        """ Sub backslashes with forwards for shell syntax"""
        path = ''
        for a in args:
                path = os.path.join(path, a)
        path = re.sub(r'\\', '/', path)
        return path


def write_key(file, value, pform=r'%.4f', ljustwidth=16):

    if type(value) != list: value = [value]
    for val in value:
        if val is None:
            file.write(str.ljust(' None', ljustwidth))
        else:
            try:
                val = float(val)
                if val < 0:
                    val = str.ljust(pform % (val), ljustwidth)
                else:
                    val = ' '+str.ljust(pform % (val), ljustwidth-1)
                file.write(val)
            except:
                file.write(str.ljust(str(val) +' ', ljustwidth-1))

def print_comment_out(text, width = 80):

    comment =" ====================== COMMENT "
    comment = re.sub(r'COMMENT', text, comment)
    for i in range(80 - len(comment)):
        comment += '='
    comment = '\n'+comment+'\n'
    sys.stdout.write(comment)

def make_xyz_matrix(matrix, nr_of_atoms):
        """ Function makes a 3d xyzmatrix from a t21-read array of geometry points """
        bohr2angs = 1.88972612499
        xyz_matrix = []
        for i in range(0, len(matrix)-nr_of_atoms, nr_of_atoms*3):
                coor = []
                for j in range(0, nr_of_atoms*3, 3):
                        coor.append([matrix[i+j]/bohr2angs, matrix[i+j+1]/bohr2angs, matrix[i+j+2]/bohr2angs])
                xyz_matrix.append(coor)
        return xyz_matrix


def read_xyz_file(file_name):
        """ reads the coordinates of a (multi)-xyz file and returns a large list, each point of which constitutes a geomety
        point of the irc. Each list elment has an x, y and z value for each atom, the total thus being a 3d-matrix. It retuns this
        matrix and the atom_list as defined in the input.

        file_name: The name (including path) for the outputfile to be analysed."""

        try:
                f = open(file_name, 'r')
        except:
                sys.stderr.write('file '+file_name+' could not be opened!\n\n')
                sys.exit()
        f_lines = f.readlines()

        conv_count = 0
        for l in f_lines:
            if re.compile(r'Geometry.*Converged', re.IGNORECASE).search(l):
                conv_count += 1

        xyz_matrix, xyz, atom_list = [], [], []
        xyz_line_found = False

        logfile = re.search(r'\.logfile', file_name)
        if conv_count < 2: logfile = False

        converged = True   ## First geometry found in IRCs is the TS-geometry.

        for line in f_lines:
            if not converged and logfile:
                converged = re.compile(r'Geometry.*Converged', re.IGNORECASE).search(line)
            elif not logfile:
                converged = True
            match = re.findall(r'\w{1,2}\b\s+.*\d+.*\s+.*\d+.*\s+.*\d.*\s', line)
            if match != []:
                temp_xyz = str.split(match[0])[1:4]
                try:
                    for j in range(3):
                        temp_xyz[j] = float(temp_xyz[j])
                    xyz.append(temp_xyz)
                    if len(xyz_matrix) == 0 and converged: atom_list.append(str.split(match[0])[0])
                    xyz_line_found = True
                except:
                    xyz_line_found = False
            else:
                xyz_line_found = False
                if xyz != [] and converged:
                    xyz_matrix.append(xyz)
                    xyz = []
                    converged = False
                elif not converged:
                    xyz = []

        if len(xyz) > 0: xyz_matrix.append(xyz)

        return [atom_list, xyz_matrix]



def new_adf(adfscript, dir, name, name_t21, atoms_block, extra, fragments = False, title=''):
        """Makes and runs a new adfscript. Takes a basic adfscript, a directory for the files, the name of the file
        and an atoms_block. the output and tape 21 files (name *.o and *.t21) are placed in the specified directory
        Optional are giving a new fragments block and extra block of input.
        adfscript: Basic adfscript read from the user input
        dir: directory were the file is placed
        name: name for the file
        name_t21: name for the tape21 of the file to be copied to dir
        atoms_block: The atoms-end block to be place in the script
        extra: Extra blocks for adfscript, this program passes an empty string if extras are not given by the user
        fragments = False: Here, a new fragment block can be passed into the adfscript for fragment analysis"""

        # finding and changing the eor statement at beginning and end of adfscript.

        new_adfscript = re.compile(r'end input', re.IGNORECASE).sub('\n'+extra, adfscript)
        for i in range(10):
            new_adfscript = re.sub(r'\n\n\n', '\n\n', new_adfscript)
        starting = """
 ====================== Starting ADF with inputscript: =========================

"""
        sys.stdout.write(starting)
        sys.stdout.write(new_adfscript)
        turbo_output = shelljoin(dir, name+'.out')
        currentPath = os.getcwd()
        dscfout  =  shelljoin(currentPath,'dscf.out')


        with open("turbo.in", 'w') as f:
            f.write(new_adfscript)
        f.close()

        with open("turbo.xyz", 'w') as f:
            for item in atoms_block:
                f.write(item)
        f.close()

        os.system('x2t turbo.xyz > coord')
        os.system('tmole')


        shutil.move(dscfout, turbo_output)
        os.remove('control')
        os.remove('statistics')
        os.remove('mos')
        os.remove('energy')
        os.remove('dscf.stat')
        os.remove('define.in')
        os.remove('define.out')
        os.remove('basis')
        os.remove('coord')
        os.remove('turbo.in')
        os.remove('turbo.xyz')





def get_input_block(specs, match_string, match_string2, default=''):
        """ reads input specifications from the input script, and returns them in a specified manner, according to
        wheter it should be string, list or two dimensional list. Defined through the type of els_return.
        specs: input specification string
        match_string: string to be matched in this case
        default: What to be returned if nothing matches. Alse determines if a list or string is returned
        match_string2: If given, a block is read and returned stripped of the match_strings
        """
        spec = default
        if re.compile(match_string, re.IGNORECASE).search(specs):
            if match_string2 and re.compile(match_string, re.IGNORECASE).search(specs):
                spec = re.compile(match_string+r'.*?'+match_string2, re.IGNORECASE|re.DOTALL).findall(specs)[0]
                spec = re.compile(match_string2+r'.*', re.IGNORECASE).sub('', spec)
                spec = re.compile(match_string+r'.*', re.IGNORECASE).sub('', spec)

        return spec

def get_input_line(specs, match_string, default=None, dim=0):
        """ reads input specifications from the input script, and returns them in a specified manner, according to
        wheter it should be string, list or two dimensional list. Defined through the type of els_return.
        specs: input specification string
        match_string: string to be matched in this case
        nr: slicing parameter for the string after splitting.
        default: What to be returned if nothing matches. Alse determines if a list or string is returned
        """
        spec = default

        p = re.compile(match_string, re.IGNORECASE)
        spec = p.findall(specs)

        if spec == []: return default

        for n in range(len(spec)):
            spec[n] = re.split(r'[ =,]*', spec[n])
            spec[n] = [elem for elem in spec[n] if elem != '']

        if dim == 0:
            spec = spec[0][0]
        elif dim == 1:
            spec = spec[0]

        return spec

class InputPrm:
    def __init__(self, name, specs='add', match_string='add', default=None, dim=0, print_names=False, unit=''):

        self.name = name
        self.unit = unit
        self.args = get_input_line(specs, match_string, default, dim)
        self.cv = None
        self.values = []
        self.columns = []

        if print_names:
            self.print_names = print_names
        elif dim and self.args is not None:
            self.print_names = []
            for arg in self.args:
                temp_name = name
                for a in arg:
                    temp_name = temp_name + '_' +a
                self.print_names.append(temp_name)
        else:
            self.print_names = [name]

        self.pwidth = 16  ## default columnwidth value!

        for name in self.print_names:
            if len(name) > self.pwidth - 2:
                self.pwidth = len(name) + 2

        self.dim = dim

    def get(self):
        return self

    def set(self, value):
        self.value = value

    def tell(self):
        """Print details of parameter"""
        print('Details of key with name "%s"' % self.name)
        print("nr of args:", self.args)
        print("names printed in data file", self.print_names)
        print("Dimension?:", self.dim)

class Parameters:
    def __init__(self):
        self.dict = {}
        self.seq = []
        self.keys_list = []

    def add(self, name, specs='Default', match_string='Default', default=None, dim=0, print_names=False, unit=''):
        self.dict[name] = InputPrm(name, specs, match_string, default, dim, print_names, unit)
        self.seq.append(self.dict[name])
        self.keys_list.append(name)
        self.__dict__[name] = self.dict[name].get()

    def __getitem__(self, name):
        return self.dict[name]

    def __setitem__(self, name, value):
        self.dict[name].set()

    def keys(self):
        return list(self.dict.keys())

    def list_keys(self):
        return self.keys_list


#-----------------------------------------------------#
#               M A I N  P R O G R A M                #
#-----------------------------------------------------#

import os, string, re, math, sys, tarfile, shutil, tempfile, time


# Description / Header
header = """
 ===============================================================================
 *  Pyfrag 2019.02                                                             *
 *  Streamlining your reaction path analysis!                                  *
 *                                                                             *
 *  Author: Willem-Jan van Zeist                                               *
 *
 *                                                                             *
 *  Find the manual at http://www.few.vu.nl/~wolters/pyfrag/                   *
 *  For some examples on how to use PyFrag, see the examples directory.        *
 *                                                                             *
 *  E-mail for PyFrag: LP.Wolters@vu.nl                                        *
 ===============================================================================
"""
sys.stdout.write(header)
print_comment_out('PyFrag Is Initializing')

## input should be: pyfrag.py jobfile [tempdir]
input = open(sys.argv[1], 'r')
adfscript = input.read()
input.close()

#seperate input specifications and adf script
specs = re.compile(r'PYFRAG\s*INPUT.*END\s*PYFRAG\s*INPUT', re.DOTALL).findall(adfscript)
if not specs: specs = re.compile(r'INPUT_?\s*SPECS.*END INPUT_SPECS', re.DOTALL).findall(adfscript)
if specs != []:
    specs = specs[0]
else:
    sys.stderr.write('Pyfrag input statments not found! Exiting Pyfrag.\n')
    sys.exit()

adfscript = re.compile(r'.*END\s*INPUT_?\s*SPECS|.*END\s*PYFRAG\s*INPUT', re.DOTALL).sub('', adfscript, 1)




# getting inline
inline = get_input_line(specs, r'inline\s(.*)', default=False)
if inline:
    try:
        inline = open(inline, 'r')
        specs += inline.read()
    except:
        sys.stderr.write('Inline file '+inline+' could not be opened.\n')

dir = os.path.dirname(sys.argv[1])
if not dir: dir = '.'

if len(sys.argv) > 2:
    temp_dir = tempfile.mkdtemp(prefix = 'Pyfrag_TMPDIR_',  dir = sys.argv[2])
else:
    temp_dir = tempfile.mkdtemp(prefix = 'Pyfrag_TMPDIR_')



# Reading wether fragment directory and/org datafile is specified
frag_dir = os.path.basename(get_input_line(specs, r'fragments?.*dir.(.*)', default='fragmentfiles'))
data_file = os.path.basename(get_input_line(specs, r'data.*file.(.*)', default='fragment_energies.txt'))

# check or make subdirectory
main_dir = dir                            # this is the main directory for the energy_file
frag_dir = shelljoin(dir,frag_dir)     # this is the dir for the fragment files

if not os.path.exists(frag_dir): os.mkdir(frag_dir)

tape21_temp_dir = shelljoin(temp_dir,'TAPE21s')

try:
    os.mkdir(tape21_temp_dir)
except:
    sys.stderr.write("Temporary TAPE21 dir is already present:"+tape21_temp_dir+"\n")

# Reading the user specs of inputfile
lt = get_input_line(specs, r'type\s*=\s*(LT)', default=False)
irc = get_input_line(specs, r'type\s*=\s*(IRC)', default=False)

if lt == False and irc == False:
    sys.stderr.write('The type of pyfrag calculation (IRC/LT/SP) has not been specified.\n')

# Reading all statements which indicate the printing of values --> printing parameters pp
pp = Parameters()

pp.add('Filename', specs,r'fa1?_name\s*=?(.*)', default='fa', unit='.out')
pp.Filename.pwidth = len(pp.Filename.args) + 7 # '_' and four numbers are added after the filenname, increase column width (plus 2 default extra)
pp.add('Point', default=[], unit='Nr.')
pp.add("Bondlength", specs, r'print.*bond\s*=?\s*(\d\s*\d.*)', dim = 2, unit='Angstrom')
pp.add("TotalIntEn", unit='kcal/mol')
pp.add("StrainFrag1", specs, r'print\s*strain\s*frag\s*1(.*)', unit='kcal/mol')
pp.add("StrainFrag2", specs, r'print\s*strain\s*frag\s*2(.*)', unit='kcal/mol')
if pp.StrainFrag1.args is not None and pp.StrainFrag2.args is not None:
    pp.add("StrainTotal", unit='kcal/mol')
    pp.add("EnergyTotal", unit='kcal/mol')

print_comment_out('Read print statements')

# Extra, if not present, they will be set as empty string and always passed on to the adfscript
extra_frag1 = get_input_block(specs, r'EXTRA\s*frag1', r'END\s*EXTRA\s*frag1')
extra_frag2 = get_input_block(specs, r'EXTRA\s*frag2', r'END\s*EXTRA\s*frag2')
extra_fa = get_input_block(specs, r'EXTRA\s*fa1?', r'END\s*EXTRA\s*fa')

# Reads specified names for fragments and fa-files
name_frag1 = get_input_line(specs, r'\s*frag1\s*=(.*)', default='frag1')
name_frag2 = get_input_line(specs, r'\s*frag2\s*=(.*)', default='frag2')

name_fa = get_input_line(specs, r'fa1?_name\s*=(.*)',default='fa')
name_fa_t21 = shelljoin(tape21_temp_dir,name_fa+'.t21')


# Defining the list of xyz_matrices for all inputs, if the frag_lists are not used (as in sp) they will be tested as false
matrix, matrix_list, matrix_list_frag1, matrix_list_frag2 = [],[],[],[]

# Read name for LT/IRC outputfile, list for the fragment atoms
# calls function to read all coordinates from the file
if lt or irc:
        output_file = get_input_line(specs, r'\s*output\s*file\s*1?\s*=?(.*)')

        output_file_2 = get_input_line(specs, r'\s*output\s*file\s*2\s*=?(.*)')

        if output_file == False and output_file_2 == False: sys.stderr.write('Outputfile(s) for irc or lt are not specified')

        LT_frag1 = get_input_block(specs, r'frag1\s*=', r'end frag1')
        LT_frag1 = re.findall(r'\b\d+\b', LT_frag1)
        for i in range(len(LT_frag1)): LT_frag1[i] = int(LT_frag1[i])

        LT_frag2 = get_input_block(specs, r'frag2\s*=', r'end frag2')
        LT_frag2 = re.findall(r'\b\d+\b', LT_frag2)
        for i in range(len(LT_frag2)): LT_frag2[i] = int(LT_frag2[i])

        #retrieving atom_list from input
        atom_list2 = get_input_block(specs, r'frag1\s*=', r'end frag2')
        atom_list2 = re.findall(r'\b\d+\b\s*\b[a-zA-Z]{1,2}\b', atom_list2)

        for i in range(len(atom_list2)):
                atom_list2[i] = re.split(' ', atom_list2[i])
                atom_list2[i][0] = int(atom_list2[i][0])
        atom_list2.sort()
        for i in range(len(atom_list2)): atom_list2[i] = atom_list2[i][1]


        if re.search(r'.xyz$|.logfile$|.amv$', output_file):
                atom_list, matrix = read_xyz_file(os.path.join(main_dir,output_file))
                if output_file_2:
                    atom_list, matrix2 = read_xyz_file(os.path.join(main_dir,output_file_2))
                    matrix.reverse()
                    matrix.extend(matrix2)


        number_of_loops = len(matrix)

        alist_length = len(atom_list)
        if len(atom_list2) < alist_length: alist_length = len(atom_list2)
        for i in range(alist_length):
            if not re.compile(atom_list2[i], re.IGNORECASE).search(atom_list[i]):
                sys.stderr.write('\nyour range of atoms does not match that of the given output file!\n')
                sys.stderr.write('Your atoms:\n')
                for j in range(len(atom_list2)):
                    sys.stderr.write(atom_list2[j]+'\n')
                sys.stderr.write('Atoms from outputfile:\n')
                for j in range(len(atom_list)):
                    sys.stderr.write(atom_list[j]+'\n')
                sys.stderr.write('\nWill continue with the atoms with atomnumbers as in the inputfile!\n')

        for xyz in matrix:
                fa_block, frag1_block, frag2_block = '', '', ''
                for j in range(len(xyz)):
                        if not j+1 in LT_frag1 and not j+1 in LT_frag2: continue
                        fa_block += atom_list[j]  + '  '
                        for k in range(len(xyz[j])):
                                fa_block += str(xyz[j][k]) + '  '
                        if j+1 in LT_frag1: fa_block += ' \n'
                        if j+1 in LT_frag2: fa_block += ' \n'
                fa_block += '\n'
                matrix_list.append(fa_block)
                for nr in LT_frag1:
                        frag1_block += atom_list[nr-1] + '  '
                        for k in range(len(xyz[nr-1])):
                                frag1_block += str(xyz[nr-1][k]) + ' '
                        frag1_block += '\n'
                frag1_block += '\n'
                matrix_list_frag1.append(frag1_block)
                for nr in LT_frag2:
                        frag2_block += atom_list[nr-1] + '  '
                        for k in range(len(xyz[nr-1])):
                                frag2_block += str(xyz[nr-1][k]) + ' '
                        frag2_block += '\n'
                frag2_block += '\n'
                matrix_list_frag2.append(frag2_block)



# For both SP and IRC/LT: get the bondlengths, angles, etc.
if pp.Bondlength.args is not None:
    for i in range(len(pp.Bondlength.args)):
        bondlengths = []
        if len(pp.Bondlength.args[i])>2:
            diffb = float(pp.Bondlength.args[i][-1])
        else:
            diffb = 0
        for xyz in matrix:
            bondlengths.append(bondlength_xyz(xyz, pp.Bondlength.args[i][0], pp.Bondlength.args[i][1]) - diffb)
        pp.Bondlength.args[i] = bondlengths


calc_points = []


data_file_initialized = False


if calc_points == []:
    for i in range(0, number_of_loops): calc_points.append(i)

pp.Point.args = []
for i in range(number_of_loops):
    pp.Point.args.append(i)
pp.Point.pwidth = len("Point") + 2
if pp.Point.pwidth < len(str(number_of_loops)) + 2: pp.Point.pwidth = len(str(number_of_loops)) + 2

print_comment_out('PyFrag Initilization done')
print_comment_out('Starting to run the ADF jobs')

# giant loop for both sp and irc/lt, for the latter it makes the fragment inputfiles
for i in calc_points:
        print_comment_out('Starting cycle: '+str(i))
        if matrix_list_frag1:
                new_name_frag1 = name_frag1+'_'+(4-len(str(i)))*'0'+str(i)
                new_name_frag1_t21 = shelljoin(frag_dir,new_name_frag1+'.out')
                new_name_frag1_t21_fa = os.path.basename(new_name_frag1+'.out')
                new_adf(adfscript, frag_dir, new_name_frag1, new_name_frag1_t21, matrix_list_frag1[i], extra_frag1, title='frag1')
#                shutil.copy(new_name_frag1_t21, '.')
        else:
                new_name_frag1_t21 = name_frag1_t21_fa
                new_name_frag1_t21_fa = name_frag1_t21_fa

        if matrix_list_frag2:
                new_name_frag2 = name_frag2+'_'+(4-len(str(i)))*'0'+str(i)
                new_name_frag2_t21 = shelljoin(frag_dir,new_name_frag2+'.out')
                new_name_frag2_t21_fa = os.path.basename(new_name_frag2+'.out')
                new_adf(adfscript, frag_dir, new_name_frag2, new_name_frag2_t21, matrix_list_frag2[i], extra_frag2, title='frag2')
#                shutil.copy(new_name_frag2_t21, '.')
        else:
                new_name_frag2_t21 = name_frag2_t21_fa
                new_name_frag2_t21_fa = name_frag2_t21_fa

        if pp.StrainFrag1.args is not None:
                pp.StrainFrag1.cv = read_energy_t21(new_name_frag1_t21) - float(pp.StrainFrag1.args)
        else:
                strain_frag1_new = False

        if pp.StrainFrag2.args is not None:
                pp.StrainFrag2.cv = read_energy_t21(new_name_frag2_t21) - float(pp.StrainFrag2.args)
        else:
                strain_frag2_new = False



        new_name_fa = name_fa+'_'+(4-len(str(i)))*'0'+str(i)
        new_name_fa_t21 = shelljoin(frag_dir,new_name_fa+'.out')

        new_adf(adfscript, frag_dir, new_name_fa, new_name_fa_t21, matrix_list[i], extra_fa, title='fragment_analysis')


        # Initializing data file if applicable
        if not data_file_initialized:
            energy_file = open(os.path.join(main_dir, data_file), 'w')
            # Print headers of columns
            for key in pp.list_keys():
                if pp[key].args is not None:
                    for name in pp[key].print_names:
                        energy_file.write(str.ljust(name, pp[key].pwidth))
            energy_file.write('\n')
            for key in pp.list_keys():
                if pp[key].args is not None:
                    for name in pp[key].print_names:
                        energy_file.write(str.ljust(pp[key].unit, pp[key].pwidth))
            energy_file.write('\n')
            energy_file.close()

        data_file_initialized = True

        pp.Filename.cv = new_name_fa

        if pp.Bondlength.args is not None:
            pp.Bondlength.cv = []
            for bond in range(len(pp.Bondlength.args)):
                pp.Bondlength.cv.append(pp.Bondlength.args[bond][i])

        pp.Point.cv = pp.Point.args[i]

        print_comment_out("PyFrag is gathering data to print")

        # pp = get_print_values(new_name_fa_t21, pp)

        if pp.StrainFrag1.args is not None and pp.StrainFrag2.args is not None:
            pp.StrainTotal.cv = pp.StrainFrag1.cv + pp.StrainFrag2.cv
            pp.EnergyTotal.cv = read_energy_t21(new_name_fa_t21) - float(pp.StrainFrag1.args) - float(pp.StrainFrag2.args)
            pp.TotalIntEn.cv  = pp.EnergyTotal.cv - pp.StrainTotal.cv
        energy_file = open(os.path.join(main_dir, data_file), 'a')
        for key in pp.list_keys():
            if pp[key].args is None: continue
            if key == "Point":
                write_key(energy_file, pp[key].cv, pform='%.0f', ljustwidth = pp[key].pwidth)
                continue
            write_key(energy_file, pp[key].cv, ljustwidth = pp[key].pwidth)
            pp[key].values.append(pp[key].cv)
        energy_file.write('\n')

        energy_file.close()
        # os.remove(os.path.join(frag_dir,'logfile'))

        print_comment_out("PyFrag finished writing to data file")


print_comment_out('PyFrag finished. Have a nice day!')

