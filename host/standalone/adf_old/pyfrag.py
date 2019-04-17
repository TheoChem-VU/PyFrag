#!/usr/bin/env python

"""
Pyfrag 2.0
Author: Willem-Jan van Zeist, vzeist@few.vu.nl.
Current Development: Lando P. Wolters, LP.Wolters@vu.nl

see also http://www.few.vu.nl/~wolters/pyfrag/ for the user manual.

This program has two main functions:
1: Reads a Linear Transit or IRC output file (normal or t21) and user defined fragments
   thereof. For each point, it generates single-point files and fragment analysis files and runs them on adf. For more information on
   the use the inputfile, check 'example_input_lt' in the examples directory. The program will generate a text file containing the decomposition energies plus other,
   user defined, values such as the strain energy.

2: Generates a number of fragment analysis files based on a single point calculation in which one variable (BONDLENGTH) is substituted
   according to a user given list of bondlenghts. The fragments are user defined by a special pyfrag input-
   file (see 'example_input_sp' in the examples directory) and the program will generate a text file containing the decomposition energies plus other, user defined,
   values.  For more information, look at the example file, it should be self-explainatory.

On a line, everything after a '#' will be ignored, allowing for comments to be inserted in your inputfile.

How to use:
Use qstc with option -pf to submit the pyfrag input script in the usual manner. Error messages
will be printed in the .err file. In the .out file, you will find the submitted adfscripts. Outputfiles of individual calculations are
stored in a subfolder (default /fragmentfiles)
After the computation, the file 'fragment_energies.txt' in the main directory will contain all the data.
"""

def get_print_values(t21, pp):
    irreps = read_irreps_t21(t21)
    frags = t21.read('Geometry', 'fragmenttype')

    pp.TotalIntEn.cv = read_energy_t21(t21)
    pp.Pauli.cv = read_pauli_t21(t21)
    pp.Elstat.cv = read_elstat_t21(t21)
    pp.Steric.cv = pp.Elstat.cv + pp.Pauli.cv
    pp.OI.cv = read_oi_t21(t21)
    if pp.IrrepOI.args is not None: pp.IrrepOI.cv = read_oi_for_irreps_t21(t21)
    if pp.PauliDecomp.args is not None: pp.PauliDecomp.cv = [read_pauli_kinetic_t21(t21), read_pauli_elstat_t21(t21)]
    if pp.VDD.args is not None: pp.VDD.cv = read_vdd_atoms_t21(t21, pp.VDD.args)
    if pp.Hirsfeld.args is not None: pp.Hirsfeld.cv = read_hirshfeld_t21(t21)
    if pp.DipoleMoment.args is not None: pp.DipoleMoment.cv = read_dipole_t21(t21)

    if pp.Pop.args is not None: pp.Pop.cv = read_populations_t21(t21, pp.Pop.args)
    if pp.OrbEn.args is not None: pp.OrbEn.cv = read_frag_orb_energies_t21(t21, pp.OrbEn.args)
    if pp.OrbEnGap.args is not None: pp.OrbEnGap.cv = read_frag_orb_energies_gap_t21(t21, pp.OrbEnGap.args)
    if pp.Overlap.args is not None: pp.Overlap.cv = read_overlaps_t21(t21, pp.Overlap.args)

    return pp

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

def read_t21(file_name1, file_name2=False, irc=True):
        """ Function which reads a tape21 (or two for irc) and returns a large list, each point of which constitutes a geomety
        point of the irc. Each list elment has an x, y and z value for each atom, the total thus being a 3d-matrix. It retuns this
        matrix and the atom_list as defined in the input.

        file_name1: The name (including path) for the tape21 to be analysed.
        file_name2 == False: A second filename for an optionally irc_backward or forward if computed seperatly, xyz_matrices will be joined together.
        IRC = True: If kept on true, IRC-tape21's will be asumed, else a LT-tape21 will be assumed.
         """
        try:
                f = kffile.kffile(file_name1)
        except:
                sys.stderr.write('file '+file_name+' could not be found!\n\n')
                sys.exit()

        nr_of_atoms = f.read('Geometry', 'nr of atoms')[0]

        #checks the input order of the atoms via a couple of lists in the t21 connecting the input and internal order of the atoms
        atomtypes = f.read('Geometry', 'atomtype') #gives the atomtypes in a list
        atomtype_index = f.read('Geometry', 'fragment and atomtype index')[-nr_of_atoms:] # designates the atom types to t21-atomnumbers
        atom_order_index = f.read('Geometry', 'atom order index')[-nr_of_atoms:] # gives nrs of atomtypes as in input
        atom_list = []
        for i in range(nr_of_atoms):
                atom_list.append([atom_order_index[i], atomtypes[atomtype_index[i]-1]])
        atom_list.sort()
        for i in range(nr_of_atoms):
                atom_list[i] = atom_list[i][1]

        if irc:
                #reading and collecting the different xyz_lists from the tape21, were they are presented with input-order atoms.
                xyz_matrix = list(f.read('IRC', 'xyz')) # This value is the TS geometry, not to be found in the IRC_forward/backward variable
                xyz_matrix = make_xyz_matrix(xyz_matrix, nr_of_atoms)

                irc_fw = f.read('IRC_Forward', 'xyz')
                irc_bw = f.read('IRC_Backward', 'xyz')
                if irc_bw is None: irc_bw = []
                if irc_fw is None: irc_fw = []
                irc_fw_mat = make_xyz_matrix(irc_fw, nr_of_atoms)
                irc_bw_mat = make_xyz_matrix(irc_bw, nr_of_atoms)
                if irc_fw != [] and irc_bw != []:
                    xyz_matrix.extend(irc_fw_mat)
                    irc_bw_mat.reverse()
                    irc_bw_mat.extend(xyz_matrix)
                    xyz_matrix = irc_bw_mat
                else:
                    xyz_matrix.extend(irc_fw_mat)
                    xyz_matrix.extend(irc_bw_mat)

                if file_name2:
                    xyz_matrix.reverse()
                    f = kffile.kffile(file_name2)
                    irc_fw = f.read('IRC_Forward', 'xyz')
                    irc_bw = f.read('IRC_Backward', 'xyz')
                    if irc_bw is None: irc_bw = []
                    if irc_fw is None: irc_fw = []
                    irc_fw_mat = make_xyz_matrix(irc_fw, nr_of_atoms)
                    irc_bw_mat = make_xyz_matrix(irc_bw, nr_of_atoms)
                    if irc_fw != [] and irc_bw != []:
                        xyz_matrix.extend(irc_fw_mat)
                        irc_bw_mat.reverse()
                        irc_bw_mat.extend(xyz_matrix)
                        xyz_matrix = irc_bw_mat
                    else:
                        xyz_matrix.extend(irc_fw_mat)
                        xyz_matrix.extend(irc_bw_mat)

        else:
                xyz_matrix = f.read('LT', 'xyz')
                xyz_matrix = make_xyz_matrix(xyz_matrix, nr_of_atoms)

        return atom_list, xyz_matrix

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

def read_gaussian_xyz_out(main_dir,output_file):
        """ reads the coordinates of an irc OR LT outputs (non-t21) and returns a file  contain all coordinate"""

        atomlist = {1:'H',3:'Li',4:'Be', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F',11:'Na',12:'Mg',13:'Al',14:'Si', 15:'P',16:'S',17:'Cl', 26:'Fe',46:'Pd'}

        file_name = shelljoin(main_dir,output_file)
        try:
                input = open(file_name, 'r')
        except:
                sys.stderr.write('file '+file_name+' could not be opened!\n\n')
                sys.exit()
        gsscript = input.read()
        input.close()

        #calculate the converged steps for every point of IRC
        step  = re.compile(r'# OF STEPS ='+r'\s*\d*', re.DOTALL).findall(gsscript)
        step  = ''.join(step)
        step  = re.compile(r'# OF STEPS =', re.DOTALL).sub('',step)
        step  = step.split()
        step  = [int(x) for x in step]
        number = 0
        loop = []
        for i in step:
           number+=i
           loop.append(number)
        #collect coordinate of each points of IRC
        data = []
        for i in loop:
            spec = []
            if re.compile(r'Input orientation:', re.IGNORECASE).search(gsscript):
                if re.compile(r'Input orientation:', re.IGNORECASE).search(gsscript):
                    testcace=re.compile(r'Input orientation:'+r'.*?'+ r'Distance matrix', re.IGNORECASE|re.DOTALL).findall(gsscript)
                    spec = re.compile(r'Input orientation:'+r'.*?'+ r'Distance matrix', re.IGNORECASE|re.DOTALL).findall(gsscript)[i]
                    spec = re.compile(r'Distance matrix (angstroms):'+r'.*', re.IGNORECASE).sub('', spec)
                    spec = re.compile(r'Input orientation:'+r'.*', re.IGNORECASE).sub('', spec)
                    newspec = []
                    spec = spec.splitlines()[5:-2]
                    atomnumber = len(spec)
                    for line in spec:
                        line = line.split()
                        line = [atomlist[int(line[1])],line[3],line[4],line[5]]
                        newspec.append(line)
                    data.append(newspec)
        #atomnumber = len(spec)
        root, ext = os.path.splitext(output_file)
        ircfile = shelljoin(main_dir,root+'.amv')
        with open(ircfile,"a") as f:
             for i in range(len(loop)):
                 f.writelines("\n")
                 for j in range(atomnumber):
                     f.writelines("\n")
                     test= data[i][j]
                     f.writelines(["%s " % item  for item in test])
        f.close()


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
        eor = re.findall(r'<<\s*[^\s>]+', adfscript)[0]
        eor = re.sub(r'[<\s]*', '', eor)
        new_adfscript = re.compile(r'ATOMS.*?END', re.DOTALL|re.IGNORECASE).sub('', adfscript)
        new_adfscript = re.sub(r'<<.*', '<<eor> '+ shelljoin(dir,name +'.out')+'\n\nATOMS\nEND', new_adfscript)
        new_adfscript = re.sub(r'[^<a-zA-Z0-9]'+eor, '\neor\n', new_adfscript)

        if re.compile(r'TITLE.*', re.IGNORECASE).search(new_adfscript):
            new_adfscript = re.compile(r'TITLE.*', re.IGNORECASE).sub('TITLE '+title, new_adfscript+'\n', 1)
        else:
            new_adfscript = re.compile(r'ATOMS', re.IGNORECASE).sub('TITLE '+title+'\n\nATOMS', new_adfscript, 1)

        if re.search('cp\s*TAPE21.*', new_adfscript):
                new_adfscript = re.compile(r'cp\s*TAPE21.*', re.IGNORECASE).sub('cp TAPE21 ' + name_t21+'\ncp logfile '+dir+'\n', new_adfscript)
        else:
                new_adfscript = re.sub(r'\neor\n', '\neor\n\ncp TAPE21 ' + name_t21+'\ncp logfile '+dir+'\n', new_adfscript)

        if re.compile(r'ATOMS.*?END', re.DOTALL|re.IGNORECASE).search(new_adfscript):
                new_adfscript = re.compile(r'ATOMS.*?END', re.DOTALL|re.IGNORECASE).sub(atoms_block + '\n'+extra, new_adfscript)
        else:
                new_adfscript = re.compile(r'end input', re.IGNORECASE).sub(atoms_block + '\n'+extra+'\nEND INPUT\n', new_adfscript)

        if fragments:
                if re.compile(r'FRAGMENTS.*?END', re.DOTALL|re.IGNORECASE).search(new_adfscript):
                        new_adfscript = re.compile(r'FRAGMENTS.*?END', re.DOTALL|re.IGNORECASE).sub(fragments, new_adfscript)
                else:
                        new_adfscript = re.compile(r'end\s*input', re.IGNORECASE).sub(fragments+'\nEND INPUT\n', new_adfscript)

        new_adfscript += '\n\nrm TAPE*\nrm logfile\n\n'

        for i in range(10):
            new_adfscript = re.sub(r'\n\n\n', '\n\n', new_adfscript)
        starting = """
 ====================== Starting ADF with inputscript: =========================

"""
        sys.stdout.write(starting)
        sys.stdout.write(new_adfscript)

        f = open('temp_pyfrag_shell_script.sh', 'w')
        f.write(new_adfscript)
        f.close()

        os.system('sh temp_pyfrag_shell_script.sh')

        os.remove('temp_pyfrag_shell_script.sh')

def make_gnuplot_script(data_file, main_dir, pp, filename):

    gnuscript = """
        # This is a gnuplot script generated by pyfrag, just run it with 'gnuplot pyfrag_gnuscript' on your command line
        # For more options of gnuplot see the website: http://www.gnuplot.info/documentation.html

        set output "pyfrag_graphs.ps"
        set terminal postscript enhanced "Times-Roman" 22
        """

    inten = str(pp.TotalIntEn.columns[0])
    pauli = str(pp.Pauli.columns[0])
    elstat = str(pp.Elstat.columns[0])
    steric = str(pp.Steric.columns[0])
    oi = str(pp.OI.columns[0])

    if pp.Bondlength.args is not None:
        x = str(pp.Bondlength.columns[0])
        xlabel = "Bondlength (\\305)"
    elif pp.Angle.args is not None:
        x = str(pp.Angle.columns[0])
        xlabel = "Angle (degrees)"
    elif pp.PlaneAngle.args is not None:
        x = str(pp.PlaneAngle.columns[0])
        xlabel = "Plane Angle (degrees)"
    elif pp.GeomVar.args is not None:
        x = str(pp.GeomVar.columns[0])
        xlabel = re.sub('_', ' ', pp.GeomVar.print_names[0])
    else:
        x = str(pp.Point.columns[0])
        xlabel = "Points"

    skip_columns = [1, 2, int(x), int(inten), int(pauli), int(elstat), int(steric), int(oi)]  # making list for the columns not automatically processed in the end of the script

    if pp.StrainFrag1.args is not None and pp.StrainFrag2.args is not None:
        strain_tot = str(pp.StrainTotal.columns[0])
        energy_tot = str(pp.EnergyTotal.columns[0])
        strain_frag1 = str(pp.StrainFrag1.columns[0])
        strain_frag2 = str(pp.StrainFrag2.columns[0])
        strain_int = [x+':'+energy_tot , x+':'+strain_tot, x+':'+inten]
        strain_decomp = [x+':'+strain_tot, x+':'+strain_frag1, x+':'+strain_frag2]
        skip_columns.extend([int(strain_tot), int(energy_tot), int(strain_frag1), int(strain_frag2)])

    decomp = [x+':'+inten, x+':'+pauli, x+':'+elstat, x+':'+oi]
    decomp_steric = [x+':'+inten, x+':'+steric, x+':'+oi]

    gnuscript += make_graph(shelljoin(main_dir, data_file), 'Decomposition Interaction Energy', \
        ['{/Symbol D}E_{int}', '{/Symbol D}E_{Pauli}', '{/Symbol D}V_{Dlstat}','E_{oi}'],decomp, xlabel=xlabel)
    gnuscript += make_graph(shelljoin(main_dir, data_file), 'Decomposition Interaction Energy', \
        ['{/Symbol D}E_{int}', '{/Symbol D}E_{Steric}','E_{oi}'],decomp_steric, xlabel=xlabel)

    if pp.StrainFrag1.args is not None and pp.StrainFrag2.args is not None:
        gnuscript += make_graph(shelljoin(main_dir, data_file), 'Strain and Interaction energy', \
            ['{/Symbol D}E_{Total}', '{/Symbol D}E_{Strain}', '{/Symbol D}E_{Int}'],strain_int, xlabel=xlabel)
        gnuscript += make_graph(shelljoin(main_dir, data_file), 'Decomposition Strain', \
            ['{/Symbol D}E_{Strain} (Total)', '{/Symbol D}E_{Strain} (Frag1)', '{/Symbol D}E_{Strain} (Frag2)'],strain_decomp, xlabel=xlabel)

    for key in pp.list_keys():
        if pp[key].args is not None:
            title = re.sub(r'_', ' ', key)
            series = []
            using = []
            for i in range(len(pp[key].print_names)):
                if pp[key].columns[i] in skip_columns:
                    continue
                series.append(re.sub(r'_', ' ', pp[key].print_names[i]))
                using.append(x+':'+str(pp[key].columns[i]))
            if series != []:
                gnuscript += make_graph(shelljoin(main_dir, data_file), title, series, using, xlabel=xlabel, ylabel=pp[key].unit)

    gnu_file = open(filename, 'w')
    gnu_file.write(gnuscript)
    gnu_file.close()


def make_graph(data_file, title, series, using, xlabel = "Bondlength (\\305)", ylabel=r"{/Symbol D}E (kcal/mol)"):
        gnuscript = """

        set rmargin 2
        set bmargin 1
        set lmargin 7
        set tmargin 0
        set ylabel "YLABEL"
        set encoding iso_8859_1
        set xlabel "XLABEL"
        set title "GRAPH_TITLE"
        plot """

        plot_line = '\"DATAFILE\" using title \'LABEL\' w lp lt -1 lw 0.5 point'
        point_styles = ['pt 1','pt 2', 'pt 3', 'pt 4', 'pt 5', 'pt 6', 'pt 7', 'pt 8']*10

        gnuscript = re.sub(r'GRAPH_TITLE', title, gnuscript)
        gnuscript = re.sub(r'XLABEL', xlabel, gnuscript)
        gnuscript = re.sub(r'YLABEL', ylabel, gnuscript)

        i = 0
        for label, use, pt in zip(series, using, point_styles):
                new_line = re.sub(r'DATAFILE', data_file, plot_line)
                new_line = re.sub(r'using', 'using ' + use, new_line)
                if label:
                    new_line = re.sub(r'LABEL', label, new_line)
                else:
                    new_line = re.sub(r'LABEL', '', new_line)
                new_line = re.sub(r'point', pt, new_line)
                i += 1
                if i < len(series):
                        new_line += ',\\\n'
                gnuscript += new_line

        return gnuscript

def do_densf(frag_t21_name, frag_orb, frag_name, cont_name):

        job = densf_job()
        job = re.sub(r'FRAGT21', frag_t21_name, job)
        job = re.sub(r'HOMOLUMO', cont_name, job)
        job = re.sub(r'CONTFILENAME', frag_name+'_cont_'+cont_name, job)
        job = re.sub(r'ORBIRREP', frag_orb[0], job)
        job = re.sub(r'ORBNUMBER', str(frag_orb[1]), job)

        os.system(job)

def densf_job():
    return """
$ADFBIN/densf<<eor>densf.out

INPUTFILE FRAGT21

DENSITY scf frag ortho

POTENTIAL coul scf

ORBITALS  SCF
ALL HOMOLUMO
END

GRID
 -3   -3  0.0
500 500
1.0  0.0   0.0 10.0
0.0  1.0   0.0 10.0
END

END INPUT
eor

$ADFBIN/cntrs<<eor>cntrs.out
SCAN
  -0.5 -0.2 -0.1 0.1 0.2 -0.5
END

DASH 0.1

FILE CONTFILENAME
SCF_ORBIRREP%ORBNUMBER 1.0

eor

rm TAPE41
"""

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

import os, string, re, math, sys, kffile, tarfile, shutil, tempfile, time

from pfm import * ## PyFragModule

# Description / Header
header = """
 ===============================================================================
 *  Pyfrag 2007.02                                                             *
 *  Streamlining your reaction path analysis!                                  *
 *                                                                             *
 *  Author: Willem-Jan van Zeist                                               *
 *  Current Development: Lando P. Wolters                                      *
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

# Removing comments from input specs.
specs = re.sub(r'#.*', '', specs)

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

try:
    os.environ['SCM_TMPDIR']
except KeyError:
    os.environ['SCM_TMPDIR'] = temp_dir

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
sp = get_input_line(specs, r'type\s*=\s*(SP)', default=False)
irc = get_input_line(specs, r'type\s*=\s*(IRC)', default=False)

if lt == False and sp == False and irc == False:
    sys.stderr.write('The type of pyfrag calculation (IRC/LT/SP) has not been specified.\n')

# Some stuff
save_tape21s = get_input_line(specs, r'save.*tape21s?')

# Reading all statements which indicate the printing of values --> printing parameters pp
pp = Parameters()

pp.add('Filename', specs,r'fa1?_name\s*=?(.*)', default='fa', unit='.out')
pp.Filename.pwidth = len(pp.Filename.args) + 7 # '_' and four numbers are added after the filenname, increase column width (plus 2 default extra)
pp.add('Point', default=[], unit='Nr.')

pp.add("Bondlength", specs, r'print.*bond\s*=?\s*(\d\s*\d.*)', dim = 2, unit='Angstrom')
pp.add("Angle", specs, r'print\s*_?angle\s*=?\s*(\d\s*\d\s*\d.*)', dim = 2, unit='Degrees')
pp.add("PlaneAngle", specs, r'print.*plane.*angle\s*=?\s*(\d\s*\d\s*\d\s*\d\s*\d.*)', dim = 2, unit='Degrees')
pp.add("AverageDistance", specs, r'print.*av.*dist.*', unit='Angstrom')
pp.add("DihedralAngle", specs, r'print\s*_?dihedralangle\s*=?\s*(\d\s*\d\s*\d\s*\d.*)', dim = 2, unit='Degrees')
pp.add("ProjectAngle", specs, r'print\s*_?projectangle\s*=?\s*(\d\s*\d\s*\d\s*\d\s*\d.*)', dim = 2, unit='Degrees')

if sp:
    pp.add("GeomVar", specs, r'geom\s*var(.*)', dim = 2, unit='Variable')
    if pp.GeomVar.args is None: pp.add("GeomVar", specs, r'Bondlengths(.*)', dim = 2, unit='Variable')

pp.add("TotalIntEn", unit='kcal/mol')
pp.add("Pauli", unit='kcal/mol')
pp.add("Elstat", unit='kcal/mol')
pp.add("Steric", unit='kcal/mol')
pp.add("OI", unit='kcal/mol')

pp.add("IrrepOI", specs, r'print\s*irreps?\s*(oi)', unit='kcal/mol')
pp.add("PauliDecomp", specs, r'print\s*\pauli\s*(decomp).*', print_names=["Pauli_Kinetic","Pauli_Elstat"], unit='kcal/mol')

pp.add("StrainFrag1", specs, r'print\s*strain\s*frag\s*1(.*)', unit='kcal/mol')
pp.add("StrainFrag2", specs, r'print\s*strain\s*frag\s*2(.*)', unit='kcal/mol')
if pp.StrainFrag1.args is not None and pp.StrainFrag2.args is not None:
    pp.add("StrainTotal", unit='kcal/mol')
    pp.add("EnergyTotal", unit='kcal/mol')

pp.add("Pop", specs, r'print.*population(.*)', dim=2, unit='Electrons')
pp.add("OrbEn", specs, r'print\s*orb.*\sen.*\s(frag.*)', dim=2, unit='eV')
pp.add("OrbEnGap", specs, r'print\s*orb.*gap(.*)', dim=2, unit='eV')
pp.add("Overlap", specs, r'print.*overlap(.*)', dim=2, unit='S')

pp.add("VDD", specs, r'print\s*VDD(.*)', dim=1, unit= 'a.u.')
pp.add("Hirsfeld", specs, r'print\s*(hirshfeld)', print_names=["HirsfeldFrag1", "HirsfeldFrag2"], unit='a.u.')
pp.add("DipoleMoment", specs, r'print\s*(dipole)', unit='Debye')

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

# Make a gnuplot script or not?
gnuplot = get_input_line(specs, r'(gnuplot)')

# Do densf?
densf = get_input_line(specs,r'densf(.*)', dim=2)

# print energies relative to point..
print_relative_energies = get_input_line(specs, r'relative.*energies.*(\d+)')

# Defining the list of xyz_matrices for all inputs, if the frag_lists are not used (as in sp) they will be tested as false
matrix, matrix_list, matrix_list_frag1, matrix_list_frag2 = [],[],[],[]

# Checking for bondlengths for sp
# Reading the fragment and total matrix
if sp:
        if pp.GeomVar.args is None:
            sys.stderr.write('\nRequired data for bondlengths is not specified!\n')
        else:
            bondlengths = pp.GeomVar.args

        bondlength_vars = []
        for i in range(len(bondlengths)):
                try:
                        bondlengths[i][0] = float(bondlengths[i][0])
                        bondlength_vars.append('BONDLENGTH')
                        number_of_loops = len(bondlengths[i])

                except ValueError:

                        if re.compile(r'step', re.IGNORECASE).search(bondlengths[i][0]) or re.compile(r'step', re.IGNORECASE).search(bondlengths[i][1]):
                                if re.compile(r'step', re.IGNORECASE).search(bondlengths[i][0]): bondlengths[i][0] = 'steps'
                                if re.compile(r'step', re.IGNORECASE).search(bondlengths[i][1]): bondlengths[i][1] = 'steps'
                                bondlengths[i].remove('steps')
                                try:
                                        bondlengths[i][0] = float(bondlengths[i][0])
                                        bondlength_vars.append('BONDLENGTH')

                                except ValueError:
                                        bondlength_vars.append(bondlengths[i][0])
                                        bondlengths[i].remove(bondlengths[i][0])

                                temp_bondlengths = []
                                try:
                                        bondlengths[i][0], bondlengths[i][1], bondlengths[i][2] = float(bondlengths[i][0]), float(bondlengths[i][1]), int(bondlengths[i][2])
                                except:
                                        sys.stderr.write('Input for bondlength steps in error, please check your input file\n')
                                        sys.stderr.write(str(bondlengths[i][0])+'  '+str(bondlengths[i][1])+'  '+str(bondlengths[i][2]))

                                if bondlengths[i][2] == 1:
                                    step = 0
                                else:
                                    step = (bondlengths[i][1] - bondlengths[i][0])/(bondlengths[i][2]-1)

                                for j in range(bondlengths[i][2]):
                                        temp_bondlengths.append(str(bondlengths[i][0]))
                                        bondlengths[i][0] += step

                                bondlengths[i] = temp_bondlengths
                        else:
                                try:
                                        bondlengths[i][0] = float(bondlengths[i][0])
                                        bondlength_vars.append('BONDLENGTH')

                                except ValueError:
                                        bondlength_vars.append(bondlengths[i][0])
                                        bondlengths[i].remove(bondlengths[i][0])

        if len(bondlengths) == 2:
            new_bondlengths = [[],[]]
            for j in range(len(bondlengths[1])):
                new_bondlengths[0].extend(bondlengths[0])
                for k in range(len(bondlengths[0])):
                    new_bondlengths[1].append(bondlengths[1][j])

            bondlengths = new_bondlengths

        number_of_loops = len(bondlengths[0])

        pp.GeomVar.args = bondlengths

        # read and makes the fragment atom_block
        fa_atoms = get_input_block(specs, r'matrix.*', r'end matrix')

        if not fa_atoms: sys.stderr.write('\nRequired matrix block for the fa-calculation is not present\n')

        for i in range(len(bondlengths[0])):
                matrix_list.append(fa_atoms)
        matrix_list = subs_geomvar(bondlength_vars, bondlengths, matrix_list)

        if re.compile(r'give\s*tape21\s*frag1\s*=?.*', re.IGNORECASE).search(specs):
                name_frag1_t21_fa = shelljoin(main_dir,get_input_line(specs, r'give\s*tape21\s*frag1\s*=?(.*)'))

        elif re.compile(r'frag1\s*=', re.IGNORECASE).search(specs):
                name_frag1_t21 = shelljoin(tape21_temp_dir,name_frag1+'.t21')
                name_frag1_t21_fa = name_frag1+'.t21'
                frag1 = get_input_block(specs, r'frag1\s*=', r'end frag1')
                if frag1 == False:
                        sys.stderr.write("frag1 is not recognized, please check your input file")

                for var in bondlength_vars:
                        if re.search(var, frag1):
                                for i in range(len(bondlengths[0])):
                                        matrix_list_frag1.append(frag1)
                                matrix_list_frag1 = subs_geomvar(bondlength_vars, bondlengths, matrix_list_frag1)
                                break

                if matrix_list_frag1 == []:
                        new_adf(adfscript, frag_dir, name_frag1, name_frag1_t21, frag1, extra_frag1, title='frag1')

        else:
                frag1_block = re.findall(r'.*f=frag1', fa_atoms)
                for i in range(1, len(frag1_block)):
                        frag1_block[0] += '\n' + frag1_block[i]
                frag1_block = 'ATOMS\n'+frag1_block[0]+'\nEND'
                frag1_block = re.sub(r'f=frag1', '', frag1_block)
                for i in range(len(bondlengths[0])):
                        matrix_list_frag1.append(frag1_block)
                matrix_list_frag1 = subs_geomvar(bondlength_vars, bondlengths, matrix_list_frag1)

        if re.compile(r'give\s*tape21\s*frag2\s*=?.*', re.IGNORECASE).search(specs):
                name_frag2_t21_fa = shelljoin(main_dir,get_input_line(specs, r'give\s*tape21\s*frag2\s*=?(.*)'))
        elif re.compile(r'frag2\s*=', re.IGNORECASE).search(specs):
                name_frag2_t21 = shelljoin(tape21_temp_dir,name_frag2+'.t21')
                name_frag2_t21_fa = name_frag2+'.t21'
                frag2 = get_input_block(specs, r'frag2\s*=', r'end frag2')
                if frag2 == False:
                        sys.stderr.write("frag2 is not recognized, please check your input file")

                for var in bondlength_vars:
                        if re.search(var, frag2):
                                for i in range(len(bondlengths[0])):
                                        matrix_list_frag2.append(frag2)
                                matrix_list_frag2 = subs_geomvar(bondlength_vars, bondlengths, matrix_list_frag2)
                                break

                if matrix_list_frag2 == []:
                        new_adf(adfscript, frag_dir, name_frag2, name_frag2_t21, frag2, extra_frag2, title='frag2')

        else:
                frag2_block = re.findall(r'.*f=frag2', fa_atoms)
                for i in range(1, len(frag2_block)):
                        frag2_block[0] += '\n' + frag2_block[i]
                frag2_block = 'ATOMS\n'+frag2_block[0]+'\nEND'
                frag2_block = re.sub(r'f=frag2', '', frag2_block)
                for i in range(len(bondlengths[0])):
                        matrix_list_frag2.append(frag2_block)
                matrix_list_frag2 = subs_geomvar(bondlength_vars, bondlengths, matrix_list_frag2)

        if pp.StrainFrag1.args is not None: strain_frag1_new = False
        if pp.StrainFrag2.args is not None: strain_frag2_new = False

        for xyz in matrix_list:
                new_xyz = []
                xyz = xyz.splitlines()[2:-1]
                for line in xyz:
                        line = line.split()
                        line = [float(line[1]), float(line[2]), float(line[3])]
                        new_xyz.append(line)
                matrix.append(new_xyz)


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

        #atom_list and xyz_matrix from outputfile
        if re.search(r'.t21$', output_file):
                if irc and not output_file_2: atom_list, matrix = read_t21(os.path.join(main_dir,output_file))
                if irc and output_file_2:
                    atom_list, matrix = read_t21(os.path.join(main_dir,output_file), os.path.join(main_dir,output_file_2))
                if lt: atom_list, matrix = read_t21(os.path.join(main_dir,output_file), irc = False)
        elif re.search(r'.xyz$|.logfile$|.amv$', output_file):
                atom_list, matrix = read_xyz_file(os.path.join(main_dir,output_file))
                if output_file_2:
                    atom_list, matrix2 = read_xyz_file(os.path.join(main_dir,output_file_2))
                    matrix.reverse()
                    matrix.extend(matrix2)
        # gaussian irc output file
        elif re.search(r'.log$|.out$', output_file):
                read_gaussian_xyz_out(main_dir,output_file)
                root, ext = os.path.splitext(output_file)
                atom_list, matrix = read_xyz_file(os.path.join(main_dir,root+'.amv'))
                if output_file_2:
                    read_gaussian_xyz_out(main_dir,output_file_2)
                    root, ext = os.path.splitext(output_file_2)
                    atom_list, matrix2 = read_xyz_file(os.path.join(main_dir,root+'.amv'))
                    matrix.reverse()
                    matrix.extend(matrix2)

        else:
                atom_list, matrix, energies = read_irc_xyz_out(os.path.join(main_dir,output_file))
                if output_file_2:
                    atom_list, matrix2, energies2 = read_irc_xyz_out(os.path.join(main_dir,output_file_2))
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
                fa_block, frag1_block, frag2_block = 'ATOMS\n', 'ATOMS\n', 'ATOMS\n'
                for j in range(len(xyz)):
                        if not j+1 in LT_frag1 and not j+1 in LT_frag2: continue
                        fa_block += atom_list[j]  + '  '
                        for k in range(len(xyz[j])):
                                fa_block += str(xyz[j][k]) + '  '
                        if j+1 in LT_frag1: fa_block += ' f=frag1\n'
                        if j+1 in LT_frag2: fa_block += ' f=frag2\n'
                fa_block += 'END\n'
                matrix_list.append(fa_block)
                for nr in LT_frag1:
                        frag1_block += atom_list[nr-1] + '  '
                        for k in range(len(xyz[nr-1])):
                                frag1_block += str(xyz[nr-1][k]) + ' '
                        frag1_block += '\n'
                frag1_block += 'END\n'
                matrix_list_frag1.append(frag1_block)
                for nr in LT_frag2:
                        frag2_block += atom_list[nr-1] + '  '
                        for k in range(len(xyz[nr-1])):
                                frag2_block += str(xyz[nr-1][k]) + ' '
                        frag2_block += '\n'
                frag2_block += 'END\n'
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

if pp.Angle.args is not None:
    for i in range(len(pp.Angle.args)):
        angles = []
        if len(pp.Angle.args[i])>3:
            diffa = float(pp.Angle.args[i][-1])
        else:
            diffa = 0
        for xyz in matrix:
            angles.append(angle_xyz(xyz, pp.Angle.args[i][0], pp.Angle.args[i][1], pp.Angle.args[i][2]) - diffa)
        pp.Angle.args[i] = angles

if pp.ProjectAngle.args is not None:
    for i in range(len(pp.ProjectAngle.args)):
        projectangles = []
        if len(pp.ProjectAngle.args[i])>5:
            diffa = float(pp.ProjectAngle.args[i][-1])
        else:
            diffa = 0
        for xyz in matrix:
            projectangles.append(projectangle_xyz(xyz, pp.ProjectAngle.args[i][0], pp.ProjectAngle.args[i][1], pp.ProjectAngle.args[i][2],pp.ProjectAngle.args[i][3],pp.ProjectAngle.args[i][4]) - diffa)
        pp.ProjectAngle.args[i] = projectangles

if pp.DihedralAngle.args is not None:
    for i in range(len(pp.DihedralAngle.args)):
        dihedralangles = []
        if len(pp.DihedralAngle.args[i])>4:
            diffa = float(pp.DihedralAngle.args[i][-1])
        else:
            diffa = 0
        for xyz in matrix:
            dihedralangles.append(dihedralangle_xyz(xyz, pp.DihedralAngle.args[i][0], pp.DihedralAngle.args[i][1], pp.DihedralAngle.args[i][2],pp.DihedralAngle.args[i][3]) - diffa)
        pp.DihedralAngle.args[i] = dihedralangles



if pp.PlaneAngle.args is not None:
    for i in range(len(pp.PlaneAngle.args)):
        plane_angles = []
        if len(pp.PlaneAngle.args[i])>5:
            diffa = float(pp.PlaneAngle.args[i][-1])
        else:
            diffa = 0
        for xyz in matrix:
            plane_angles.append(plane_angle_xyz(xyz, pp.PlaneAngle.args[i][0], pp.PlaneAngle.args[i][1], pp.PlaneAngle.args[i][2], pp.PlaneAngle.args[i][3], pp.PlaneAngle.args[i][4]) - diffa)
        pp.PlaneAngle.args[i] = plane_angles

if pp.AverageDistance.args is not None:
    pp.AverageDistance.args = average_distance(matrix)


starting_point = int(get_input_line(specs, r'starting.*point(.*)', default=0))
points_range = get_input_line(specs, r'points?.*range(.*)', dim=1, default=[])

calc_points = []

for entry in points_range:

    if re.compile(r'odd', re.IGNORECASE).search(entry):
        for i in range(1, number_of_loops, 2):
            if i >= starting_point: calc_points.append(i)
    elif re.compile(r'even', re.IGNORECASE).search(entry):
        for i in range(0, number_of_loops, 2):
            if i >= starting_point: calc_points.append(i)
    else:
        temp = re.split(r'-', entry)

        if temp[0] == '': temp[0] = 0
        if temp[-1] == '': temp[-1] = number_of_loops

        if len(temp)>1:
            for i in range(int(temp[0]), int(temp[-1])+1):
                if i >= starting_point: calc_points.append(i)
        else:
            if i >= starting_point: calc_points.append(int(temp[0]))

data_file_initialized = False
if starting_point != 0:
    data_file_initialized = True

if calc_points == []:
    for i in range(starting_point, number_of_loops): calc_points.append(i)

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
                new_name_frag1_t21 = shelljoin(tape21_temp_dir,new_name_frag1+'.t21')
                new_name_frag1_t21_fa = os.path.basename(new_name_frag1+'.t21')
                new_adf(adfscript, frag_dir, new_name_frag1, new_name_frag1_t21, matrix_list_frag1[i], extra_frag1, title='frag1')
                shutil.copy(new_name_frag1_t21, '.')
        else:
                new_name_frag1_t21 = name_frag1_t21_fa
                new_name_frag1_t21_fa = name_frag1_t21_fa

        if matrix_list_frag2:
                new_name_frag2 = name_frag2+'_'+(4-len(str(i)))*'0'+str(i)
                new_name_frag2_t21 = shelljoin(tape21_temp_dir,new_name_frag2+'.t21')
                new_name_frag2_t21_fa = os.path.basename(new_name_frag2+'.t21')
                new_adf(adfscript, frag_dir, new_name_frag2, new_name_frag2_t21, matrix_list_frag2[i], extra_frag2, title='frag2')
                shutil.copy(new_name_frag2_t21, '.')
        else:
                new_name_frag2_t21 = name_frag2_t21_fa
                new_name_frag2_t21_fa = name_frag2_t21_fa

        if pp.StrainFrag1.args is not None:
                frag1_t21 = kffile.kffile(new_name_frag1_t21)
                pp.StrainFrag1.cv = read_energy_t21(frag1_t21) - float(pp.StrainFrag1.args)
        else:
                strain_frag1_new = False

        if pp.StrainFrag2.args is not None:
                frag2_t21 = kffile.kffile(new_name_frag2_t21)
                pp.StrainFrag2.cv = read_energy_t21(frag2_t21) - float(pp.StrainFrag2.args)
        else:
                strain_frag2_new = False

        fragment_block = 'FRAGMENTS \nfrag1  ' + new_name_frag1_t21_fa + '\nfrag2  ' + new_name_frag2_t21_fa + '\nEND\n'

        new_name_fa = name_fa+'_'+(4-len(str(i)))*'0'+str(i)
        new_name_fa_t21 = shelljoin(tape21_temp_dir,new_name_fa+'.t21')

        new_adf(adfscript, frag_dir, new_name_fa, new_name_fa_t21, matrix_list[i], extra_fa, fragment_block, title='fragment_analysis')
        if matrix_list_frag1: os.remove(new_name_frag1_t21_fa)
        if matrix_list_frag2: os.remove(new_name_frag2_t21_fa)

        current_fa_t21 = kffile.kffile(new_name_fa_t21)

        irreps = remove_sub_irreps(read_irreps_t21(current_fa_t21))
        pp.IrrepOI.print_names = ["OI_irrep_"+irreps[0]]
        for ir in range(1, len(irreps)):
            pp.IrrepOI.print_names.append("OI_irrep_"+irreps[ir])


        # Initializing data file if applicable
        if not data_file_initialized:
            energy_file = open(os.path.join(main_dir, data_file), 'w')
            # Print headers of columns
            for key in pp.list_keys():
                if pp[key].args is not None:
                    for name in pp[key].print_names:
#                        energy_file.write(string.ljust(name, pp[key].pwidth))
                        energy_file.write(str.ljust(name, pp[key].pwidth))
            energy_file.write('\n')
            # Print second line: units
            for key in pp.list_keys():
                if pp[key].args is not None:
                    for name in pp[key].print_names:
                        energy_file.write(str.ljust(pp[key].unit, pp[key].pwidth))
#                        energy_file.write(string.ljust(pp[key].unit, pp[key].pwidth))
            energy_file.write('\n')
            energy_file.close()

        data_file_initialized = True

        pp.Filename.cv = new_name_fa

        if pp.Bondlength.args is not None:
            pp.Bondlength.cv = []
            for bond in range(len(pp.Bondlength.args)):
                pp.Bondlength.cv.append(pp.Bondlength.args[bond][i])

        if pp.Angle.args is not None:
            pp.Angle.cv = []
            for bond in range(len(pp.Angle.args)):
                pp.Angle.cv.append(pp.Angle.args[bond][i])

        if pp.DihedralAngle.args is not None:
            pp.DihedralAngle.cv = []
            for bond in range(len(pp.DihedralAngle.args)):
                pp.DihedralAngle.cv.append(pp.DihedralAngle.args[bond][i])

        if pp.ProjectAngle.args is not None:
            pp.ProjectAngle.cv = []
            for bond in range(len(pp.ProjectAngle.args)):
                pp.ProjectAngle.cv.append(pp.ProjectAngle.args[bond][i])

        if pp.PlaneAngle.args is not None:
            pp.PlaneAngle.cv = []
            for bond in range(len(pp.PlaneAngle.args)):
                pp.PlaneAngle.cv.append(pp.PlaneAngle.args[bond][i])

        if pp.AverageDistance.args is not None:
            pp.AverageDistance.cv = pp.AverageDistance.args[i]

        if sp:
            if pp.GeomVar is not None:
                pp.GeomVar.cv = []
                for gv in range(len(pp.GeomVar.args)):
                    pp.GeomVar.cv.append(pp.GeomVar.args[gv][i])

        pp.Point.cv = pp.Point.args[i]

        print_comment_out("PyFrag is gathering data to print")

        pp = get_print_values(current_fa_t21, pp)

        if pp.StrainFrag1.args is not None and pp.StrainFrag2.args is not None:
            pp.StrainTotal.cv = pp.StrainFrag1.cv + pp.StrainFrag2.cv
            pp.EnergyTotal.cv = pp.StrainFrag1.cv + pp.StrainFrag2.cv + pp.TotalIntEn.cv

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
        os.remove(os.path.join(frag_dir,'logfile'))

        print_comment_out("PyFrag finished writing to data file")
        if densf:
            print_comment_out("Creating some orbitals plots with densf and cntrs")
            for orb in densf:
                if orb[0] == 'frag1':
                    frag_orb = find_holu_t21(current_fa_t21, get_frag_nr(current_fa_t21, orb[0]), orb[1], frag_nrs=True)
                    do_densf(new_name_frag1_t21, frag_orb, shelljoin(frag_dir, new_name_frag1), orb[1])
                if orb[0] == 'frag2':
                    frag_orb = find_holu_t21(current_fa_t21, get_frag_nr(current_fa_t21, orb[0]), orb[1], frag_nrs=True)
                    do_densf(new_name_frag2_t21, frag_orb, shelljoin(frag_dir, new_name_frag2), orb[1])
                if orb[0] == 'fa':
                    fa_orb = find_holu_t21_fa(current_fa_t21, orb[1])
                    do_densf(new_name_fa_t21, fa_orb, shelljoin(frag_dir, new_name_fa), orb[1])

        if save_tape21s is None:
            if not re.compile(r'give\s*tape21\s*frag1\s*=?.*', re.IGNORECASE).search(specs):
                os.remove(new_name_frag1_t21)
            if not re.compile(r'give\s*tape21\s*frag2\s*=?.*', re.IGNORECASE).search(specs):
                os.remove(new_name_frag2_t21)
            os.remove(new_name_fa_t21)

if save_tape21s:
    tape21_save_dir = shelljoin(frag_dir, 'TAPE21s')
    try:
        os.mkdir(tape21_save_dir)
    except:
        sys.stderr.write("Dir for saving TAPE21s already present :"+tape21_save_dir+"\n")
    for file in os.listdir(tape21_temp_dir):
        shutil.copyfile(os.path.join(tape21_temp_dir, file), os.path.join(tape21_save_dir, file))

os.chdir(main_dir)
shutil.rmtree(temp_dir, ignore_errors = True)

columns = 1
column_keys = {}
for key in pp.list_keys():
    if pp[key].args is not None:
        for name in pp[key].print_names:
            pp[key].columns.append(columns)
            column_keys[columns] = key
            columns += 1

if gnuplot is not None:
    print_comment_out('Generating GnuPlot script file')
    make_gnuplot_script(data_file, main_dir, pp, sys.argv[1] + "_pyfrag_gnuscript")

if print_relative_energies is not None:
    energy_file = open(os.path.join(main_dir, data_file), 'r')
    energy_file = energy_file.readlines()
    headers, units = energy_file[0], energy_file[1]
    data = [line.split() for line in energy_file[2:]]
    #data = data[2:]
    data_cols = [[str(d[0]) for d in data]]
    data_cols.append([str(d[1]) for d in data])
    for i in range(2, columns-1):
        data_cols.append([float(d[i]) for d in data])
        if column_keys[i+1] != "GeomVar":
            data_cols[-1] = [nr - data_cols[-1][int(print_relative_energies)] for nr in data_cols[-1]]

    refile = open(os.path.join(main_dir, 'relative_'+data_file), 'w')
    refile.write(headers+units)
    for i in range(len(data_cols[0])):
        for j in range(columns - 1):
            write_key(refile, data_cols[j][i], ljustwidth = pp[column_keys[j+1]].pwidth)
        refile.write('\n')

    if gnuplot is not None:
        make_gnuplot_script('relative_'+data_file, main_dir, pp, sys.argv[1] + "_relative_pyfrag_gnuscript")


print_comment_out('PyFrag finished. Have a nice day!')
