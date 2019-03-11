###############
# PyFragModule, some modules that are used in the pyfrag program
###############

import string, os, sys, re, math

def read_energy_t21(t21, unit='kcal/mol'):

    energy = t21.read('Energy', 'Bond Energy')[0]
    if unit == 'kcal/mol': energy = 627.50954*energy

    return energy

def read_pauli_t21(t21, unit='kcal/mol'):

    energy = t21.read('Energy', 'Pauli Total')[0]
    if unit == 'kcal/mol': energy = 627.50954*energy

    return energy

def read_elstat_t21(t21, unit='kcal/mol'):

    energy = t21.read('Energy', 'Electrostatic Interaction')[0]
    if unit == 'kcal/mol': energy = 627.50954*energy

    return energy

def read_oi_t21(t21, unit='kcal/mol'):

    energy = t21.read('Energy', 'Orb.Int. Total')[0]
    if unit == 'kcal/mol': energy = 627.50954*energy

    return energy

def read_oi_for_irreps_t21(t21, unit='kcal/mol'):

    irreps = read_irreps_t21(t21)
    irreps = remove_sub_irreps(irreps)
    irreps_oi = []
    total_oi_ir = 0
    energy_convert = 627.50954

    for irrep in irreps:
        print (irrep)
        read_oi = t21.read('Energy', 'Orb.Int. '+irrep)[0]
        if read_oi or read_oi == 0.0:
            irreps_oi.append(read_oi)
            total_oi_ir += read_oi
        else:           
            sys.stderr.write('\nWrong irrep input for orbital interaction irrep.\n')
            irreps_oi.append(None)

    # scaling to for the fit-correction, which is incorporated in the total orb int, but not in the irreps thereof, in the output it is simply scaled
    fitcorr = t21.read('Energy', 'Orb.Int. Total')[0]/total_oi_ir

    for i in range(len(irreps_oi)):
        irreps_oi[i] = irreps_oi[i]*fitcorr*energy_convert

    return irreps_oi

def read_pauli_kinetic_t21(t21):
    au2kcal = 627.50954
    return au2kcal*t21.read('Energy', 'Pauli Kinetic')[0]

def read_pauli_elstat_t21(t21):
    au2kcal = 627.50954
    return au2kcal*(t21.read('Energy', 'Pauli Total')[0] - t21.read('Energy', 'Pauli Kinetic')[0])

def read_irreps_t21(t21):
    irreps = t21.read('Symmetry', 'symlab')
    return irreps

def remove_sub_irreps(irreps):
    new_irreps = []
    for i in irreps:
        irrep = re.sub(r':.*', '', i)
        if not irrep in new_irreps: new_irreps.append(irrep)
    return new_irreps

def read_irreps_frag_t21(t21, frag_nr):
    irreps = t21.read('Ftyp '+str(frag_nr), 'bb')
    return irreps

def nr_core_orbs_irrep_t21(t21, irrep):

    nocore_mo_count = t21.read(irrep, 'nmo_A')[0]
#    including_core_mo_count =len(t21.read(irrep, 'Total SFO Gross Populations in t'))
        
#    return including_core_mo_count - nocore_mo_count
    return 0   
def read_irreps_no_subspecies_t21(t21):

    irreps = t21.read('Symmetry', 'symlab')

    irreps_done = []
    for irrep in irreps:
        irrep = re.sub(r'[.:].*', '', irrep)
        if irrep in irreps_done:
            continue
        else:
            irreps_done.append(irrep)

    return irreps_done

def read_overlap_t21(t21, omat, orb1, orb2):

    if omat is None:
        sys.stderr.write("Something's wrong with overlap orbital irrep!\n")
        return None

    try:
        orb1 = int(orb1)
        orb2 = int(orb2)
    except:
        sys.stderr.write("something's wrong with overlap orbital nrs!\n")
        return None
    
    #making sure orb1 has the larger value (necessary below)
    if orb2 > orb1:    
        temp = orb2
        orb2 = orb1
        orb1 = temp
            
    # finding the right orbital number in the triangular matrix omat (1,1)(2,1)(2,2)(3,1)...
    # since this matrix is printed in one array, below it first finds the right row, which is
    # the overlap_nr at the end of the for loop, positioning it at the end of the row, the right
    # spot on that row is then defined through orb2 
    overlap_nr = 0
    for i in range(orb1):
        overlap_nr += orb1 - i 
    overlap_nr = overlap_nr - orb1 + orb2

    return abs(omat[overlap_nr-1])

def read_overlaps_t21(t21, orbs):

    ol_list = []
    omat_dict = {}
    orb_list = []
    for orb in orbs:
        if len(orb) < 3 and re.search(r'HOMO|LUMO', orb[-1]):
            orb_1 = string.split(orb[-1], '|')[0]
            orb_2 = string.split(orb[-1], '|')[1]
            print(orb_1, orb_2)
            irrep1, orb_nr1 = find_holu_t21(t21, get_frag_nr(t21, 'frag1'), orb_1)
            irrep2, orb_nr2 = find_holu_t21(t21, get_frag_nr(t21, 'frag2'), orb_2)
            
            if irrep1 == None or irrep2 == None:
                sys.stderr.write("requested HOMO/LUMO for overlap value does not seem to exist\n")
                ol_list.append(None)
                continue
            
            if irrep1 != irrep2:
                sys.stderr.write("requested HOMO/LUMO orbital overlap values --> orbs do not have same symmetry!\n")
                ol_list.append(0)
                continue
            else:
                irrep = irrep1

        else:
            print(orb[-2], orb[-1])
            irrep = orb[-2] 
            orb_nr1 = string.split(orb[-1], '|')[0]
            orb_nr2 = string.split(orb[-1], '|')[1]
#            irrep, orb_nr1, orb_nr2 = orb[-3], orb[-2], orb[-1]
#            irrep = orb[-3], orb_nr1 = orb[-2], orb_nr2 = orb[-1]
        if not irrep in omat_dict: omat_dict[irrep] = t21.read(irrep, 'S-CoreSFO')
        ol_list.append(read_overlap_t21(t21, omat_dict[irrep], orb_nr1, orb_nr2))
        
    return ol_list
    
def read_sfo_orbs_t21(t21, irreps):
    SFO_energies = []
    SFO_occ = []
    SFO_list = []
    
    #makking two lists for occupations and energies
    for sym in irreps:
        SFO_energies = list(t21.read(sym, 'eps_A'))
        SFO_occ = list(t21.read(sym, 'froc_A'))

        for i in range(len(SFO_energies)):
            SFO_list.append([SFO_energies[i], SFO_occ[i], sym, i+1, 'spin A'])

        if t21.read(sym, 'froc_B') != None:
            SFO_occ = list(t21.read(sym, 'froc_B'))
            SFO_energies = list(t21.read(sym, 'eps_B'))
            for i in range(len(SFO_energies)):
                SFO_list.append([SFO_energies[i], SFO_occ[i], sym, i+1, 'spin B'])
    
    SFO_list.sort()                          
    return SFO_list

def read_frag_orb_energy_t21(t21, frag_nr, irrep, orb_nr):

    au2ev = 27.212

    energy = t21.read('Ftyp '+str(frag_nr)+irrep, 'eps')
    if energy == None:
        # removing .u or .g since they shouldn't be included for this variable
        new_irrep = re.sub(r'\.g|\.u', '', irrep)
        energy = t21.read('Ftyp 1'+new_irrep, 'eps')
        if energy == None:
            sys.stderr.write('\nWrong irreps or orbital numbers for the orbital energies.\n')
            return None
        
    return au2ev*energy[int(orb_nr) - 1]

def read_frag_orb_energies_t21(t21, orbs):
    
    en_list = []
    for orb in orbs:
        frag_nr = get_frag_nr(t21, orb[0])
        if len(orb) < 3 and re.search(r'HOMO|LUMO', orb[-1]):
            irrep, orb_nr = find_holu_t21(t21, frag_nr, orb[-1], frag_nrs = True)
        else:
            irrep, orb_nr = orb[1], orb[2]
            
        en_list.append(read_frag_orb_energy_t21(t21, frag_nr, irrep, orb_nr))

    return en_list

def read_frag_orb_energies_gap_t21(t21, orbs):

    en_list = []
    frag1 = get_frag_nr(t21, 'frag1')
    frag2 = get_frag_nr(t21, 'frag2')
    for orb in orbs:
        if re.search(r'\|', orb[0]) and len(orb) == 1:
            orb = orb[0].split('|')
        irrep1, orb_nr1 = find_holu_t21(t21, frag1, orb[0], frag_nrs = True)
        irrep2, orb_nr2 = find_holu_t21(t21, frag2, orb[1], frag_nrs = True)
        
        e1 = read_frag_orb_energy_t21(t21, frag1, irrep1, orb_nr1)
        e2 = read_frag_orb_energy_t21(t21, frag2, irrep2, orb_nr2)

        en_list.append(e1 - e2)

    return en_list

def get_frag_nr(t21, frag):
    
    frags = t21.read('Geometry', 'fragmenttype')
        
    if frags[0] == frag:
        return 1
    elif frags[1] == frag:
        return 2

def find_holu_index(holu):
    
    if re.search(r'HOMO', holu):
        homo_index = re.sub(r'HOMO', '', holu)
        if homo_index == '': 
            homo_index = 0
        else:
            try:
                homo_index = int(homo_index)
            except:
                homo_index = 0

    if re.search(r'LUMO', holu):
        homo_index = re.sub(r'LUMO', '', holu)
        if homo_index == '':
            homo_index = 0 + 1
        else:
            try:
                homo_index = int(homo_index) + 1
            except:
                homo_index = 0 + 1
    return homo_index

def find_holu_t21_fa(t21, holu):
    homo_index = find_holu_index(holu)
    irreps = read_irreps_t21(t21)
    orbs = []
    for irrep in irreps:
        temp_orbs = []
        occs = t21.read(irrep, 'froc_A')
        ens = t21.read(irrep, 'eps_A')
        for i in range(len(occs)):
            temp_orbs.append([ens[i], occs[i], irrep, i+1])
        orbs.extend(temp_orbs)

    orbs.sort()

    for i in range(1, len(orbs)):
        if orbs[-i][1] == 0.0 and orbs[-i-1][1] != 0.0:
            holu_orb = [orbs[-i-1+homo_index][-2], orbs[-i-1+homo_index][-1]]  # Orbital irrep & nr in the fragment analysis calc.

    return holu_orb


def find_holu_t21(t21, frag_nr, holu, frag_nrs = False):

    homo_index = find_holu_index(holu)

    frag_irreps = read_irreps_frag_t21(t21, frag_nr)
    
    frag_orb_nrs = t21.read('SFOs', 'ifo')
    orb_nrs_fa = t21.read('SFOs', 'isfo')
    irreps_frags = t21.read('SFOs', 'subspecies')
    fragment_nr = t21.read('SFOs', 'fragment')
    energies = t21.read('SFOs', 'energy')
    occs = t21.read('SFOs', 'occupation')
    
    irreps = read_irreps_t21(t21)
    nr_of_SFOs_per_irrep = t21.read('Symmetry', 'norb')
    irreps_fa = []
    for i, nr in zip(irreps, nr_of_SFOs_per_irrep):
        for orb in range(1, nr+1):
            irreps_fa.append(i)
    
    frag_orbs = []
    for e, o, i, f, fn, onf, irf in zip(energies, occs, irreps_frags, frag_orb_nrs, fragment_nr, orb_nrs_fa, irreps_fa):
        if re.search('HOMO', holu) and o == 0.0: continue
        if re.search('LUMO', holu) and o != 0.0: continue
        if fn == int(frag_nr): frag_orbs.append([e, o, i, f, fn, onf, irf])
    frag_orbs.sort()

    if re.search('HOMO', holu):
        if len(frag_orbs) - 1 + homo_index < 0: return None, None
        holu_orb = [frag_orbs[-1+homo_index][-1], frag_orbs[-1+homo_index][-2]]  # Orbital irrep & nr in the fragment analysis calc.
        holu_orb_frag = [frag_orbs[-1+homo_index][-5], frag_orbs[-1+homo_index][-4]]

    if re.search('LUMO', holu):
        holu_orb = [frag_orbs[homo_index - 1][-1], frag_orbs[homo_index - 1][-2]]  # Orbital irrep & nr in the fragment analysis calc.
        holu_orb_frag = [frag_orbs[homo_index - 1][-5], frag_orbs[homo_index - 1][-4]]

    holu_orb[1] += nr_core_orbs_irrep_t21(t21, holu_orb[0])
    
    if frag_nrs:
        return holu_orb_frag[0], holu_orb_frag[1]
    else:
        return holu_orb[0], holu_orb[1]

def read_population_t21(t21, irrep, orb_nr):
    if re.search('frag', irrep) and re.search(r'HOMO|LUMO', orb_nr):
        frag_nr = get_frag_nr(t21, irrep)
        irrep, orb_nr = find_holu_t21(t21, frag_nr, orb_nr)
        print(irrep,orb_nr) 
    if irrep == None: 
        sys.stderr.write('\nWrong irreps.\n')
        return None
    #popul = t21.read(irrep, 'SFO')
    popul = t21.read('SFO popul', 'sfo_grosspop')
    if popul != None:
       return popul[int(orb_nr)-1]
       #return popul[int(orb_nr)]
    else:
        # removing .u or .g since they shouldn't be included for this variable
        new_irrep = re.sub(r'\.g|\.u', '', irrep)
        #popul = t21.read(new_irrep, 'SFO')
        popul = t21.read('SFO popul', 'sfo_grosspop')
        if popul != None and popul!=[]:
            return popul[int(orb_nr)-1]
            #return popul[int(orb_nr)]
        else:
            sys.stderr.write('\nWrong populations.\n')
            return None

def read_populations_t21(t21, pops):
    pop_list = []
    for pop in pops:
        if re.search('frag', pop[-2]) and re.search(r'5d|d5', pop[-1]):
            popul = 0
            for orb in ['LUMO','LUMO-1','LUMO-2','HOMO','HOMO-1','HOMO-2','HOMO-3','HOMO-4']:
                popul += read_population_t21(t21, pop[-2], orb)
            pop_list.append(popul)
        else:
            pop_list.append(read_population_t21(t21, pop[-2], pop[-1]))
    return pop_list
                 
def read_vdd_atoms_t21(t21, atom_list):

    vdd_list = []
    for atom_nr in atom_list:
        try:
            vdd_scf = t21.read('Properties', 'AtomCharge_SCF Voronoi')[int(atom_nr)-1]
            vdd_init = t21.read('Properties', 'AtomCharge_initial Voronoi')[int(atom_nr)-1]
            vdd_list.append(vdd_scf - vdd_init)
        except:
            sys.stderr.write('\nSomethings wrong with vdd_atoms: '+str(atom_nr)+'\n')
            vdd_list.append(None)
    return vdd_list

def read_hirshfeld_t21(t21):
    
    hirshfeld = t21.read('Properties', 'FragmentCharge Hirshfeld')
    hirshfeld_frag1 = hirshfeld[get_frag_nr(t21, 'frag1') - 1]
    hirshfeld_frag2 = hirshfeld[get_frag_nr(t21, 'frag2') - 1]
    return [hirshfeld_frag1, hirshfeld_frag2]

def read_dipole_t21(t21):

    au2debye = 2.54177
    dipole = t21.read('Properties', 'Dipole')
    dipole = au2debye*math.sqrt(dipole[0]**2 + dipole[1]**2 + dipole[2]**2)
    return dipole

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

def angle_xyz(xyz, atom1, atom2, atom3):

    try:
        angle = [int(atom1), int(atom2), int(atom3)]
    except:
        sys.stderr.write("values for printing bondlength are not integers!")
        return None

    vec1 = [xyz[angle[0]-1][0] - xyz[angle[1]-1][0], xyz[angle[0]-1][1] - xyz[angle[1]-1][1], xyz[angle[0]-1][2] - xyz[angle[1]-1][2]]
    vec2 = [xyz[angle[2]-1][0] - xyz[angle[1]-1][0], xyz[angle[2]-1][1] - xyz[angle[1]-1][1], xyz[angle[2]-1][2] - xyz[angle[1]-1][2]]

    return determine_angle(vec1, vec2)


# angle between the bond and a plane
def projectangle_xyz(xyz, atom1, atom2, atom3,atom4,atom5):
    
    try:
        angle = [int(atom1), int(atom2), int(atom3),int(atom4),int(atom5)]
    except:
        sys.stderr.write("values for printing dihedral are not integers!")
        return None
    vec1 = [xyz[angle[0]-1][0] - xyz[angle[1]-1][0], xyz[angle[0]-1][1] - xyz[angle[1]-1][1], xyz[angle[0]-1][2] - xyz[angle[1]-1][2]]
    vec2 = [xyz[angle[2]-1][0] - xyz[angle[1]-1][0], xyz[angle[2]-1][1] - xyz[angle[1]-1][1], xyz[angle[2]-1][2] - xyz[angle[1]-1][2]]
    vec3 = [xyz[angle[4]-1][0] - xyz[angle[3]-1][0], xyz[angle[4]-1][1] - xyz[angle[3]-1][1], xyz[angle[4]-1][2] - xyz[angle[3]-1][2]]
    n_vec1 = vector_outerproduct(vec1,vec2)
    m_vec1 = vector_normalize(n_vec1)
    m_vec2 = vector_normalize(vec3)
    dot = vector_innerproduct(m_vec1, m_vec2)

    angle = math.acos(dot)

    return 90-angle*180/math.pi

def dihedralangle_xyz(xyz, atom1, atom2, atom3,atom4):


    try:
        angle = [int(atom1), int(atom2), int(atom3),int(atom4)]
    except:
        sys.stderr.write("values for printing dihedral are not integers!")
        return None
    vec1 = [xyz[angle[0]-1][0] - xyz[angle[1]-1][0], xyz[angle[0]-1][1] - xyz[angle[1]-1][1], xyz[angle[0]-1][2] - xyz[angle[1]-1][2]]
    vec2 = [xyz[angle[2]-1][0] - xyz[angle[1]-1][0], xyz[angle[2]-1][1] - xyz[angle[1]-1][1], xyz[angle[2]-1][2] - xyz[angle[1]-1][2]]
    vec3 = [xyz[angle[1]-1][0] - xyz[angle[2]-1][0], xyz[angle[1]-1][1] - xyz[angle[2]-1][1], xyz[angle[1]-1][2] - xyz[angle[2]-1][2]]
    vec4 = [xyz[angle[3]-1][0] - xyz[angle[2]-1][0], xyz[angle[3]-1][1] - xyz[angle[2]-1][1], xyz[angle[3]-1][2] - xyz[angle[2]-1][2]]
    n_vec1 = vector_outerproduct(vec1,vec2)
    n_vec2 = vector_outerproduct(vec3,vec4)
    m_vec1 = vector_normalize(n_vec1)
    m_vec2 = vector_normalize(n_vec2)
    dot = vector_innerproduct(m_vec1, m_vec2)

    angle = math.acos(dot)

    return angle*180/math.pi


def plane_angle_xyz(xyz, atom1, atom2, atom3, atom4, atom5):
    
    atom1, atom2, atom3, atom4, atom5 = int(atom1), int(atom2), int(atom3), int(atom4), int(atom5)
    xyz = set_origin_xyz(xyz, atom3)
    
    n_vec1 = vector_outerproduct(xyz[atom1-1], xyz[atom2-1])
    n_vec2 = vector_outerproduct(xyz[atom4-1], xyz[atom5-1])
    
    angle = determine_angle(n_vec1, n_vec2)    
    return angle

def vector_normalize(vec):

    s = math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    if s != 0.0:
        vec[0] = vec[0]/s
        vec[1] = vec[1]/s
        vec[2] = vec[2]/s
    return vec

def vector_innerproduct(vec1, vec2):
    
    inner = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]
    return inner

def vector_outerproduct(vec1, vec2):
    x = vec1[1]*vec2[2] - vec2[1]*vec1[2]
    y = -vec1[0]*vec2[2] + vec2[0]*vec1[2]
    z = vec1[0]*vec2[1] - vec2[0]*vec1[1]
    return [x, y, z]

def determine_angle(vec1, vec2):
    
    vec1 = vector_normalize(vec1)
    vec2 = vector_normalize(vec2)
    dot = vector_innerproduct(vec1, vec2)

    angle = math.acos(dot)

    return angle*180/math.pi    

def set_origin_xyz(xyz, atom):
    import os, sys, re, math, string

    origin = int(atom) - 1
    xc = xyz[origin][0]
    yc = xyz[origin][1]
    zc = xyz[origin][2]
    for i in range(len(xyz)):

        xyz[i][0] = xyz[i][0] - xc
        xyz[i][1] = xyz[i][1] - yc
        xyz[i][2] = xyz[i][2] - zc

    return xyz

def read_irc_xyz_out(file_name):
        """ reads the coordinates of an irc OR LT outputs (non-t21) and returns a large list, each point of which constitutes a geomety
        point of the irc. Each list elment has an x, y and z value for each atom, the total thus being a 3d-matrix. It retuns this
        matrix and the atom_list as defined in the input. Also returns energies relative to starting point.

        file_name: The name for the outputfile to be analyzed."""

        import re, string, sys
        
        try:
                f = open(file_name, 'r')
        except:
                sys.stderr.write('file '+file_name+' could not be opened!\n\n')
                sys.exit()
        f = f.read()
        
        #Finding the atom list
        atoms = re.compile(r'ATOMS.*?=====\s+X Y Z.*?FRAGMENTS', re.DOTALL).findall(f)[0]
        atoms = re.findall(r'\s+\d+\s{2}\w{1,2}\s+-?\d\.\d{4}', atoms)
        for i in range(len(atoms)):
                atoms[i] = re.findall(r'[a-zA-Z]{1,2}', atoms[i])[0]

        summary = re.compile(r'Summary of the IRC.*?End of Path Summary', re.DOTALL).findall(f)
        if summary == []:
            summary = re.compile(r'Summary of the LT.*?End of Path Summary', re.DOTALL).findall(f)
        if summary == []:
            return
        else:
            summary = summary[0]

        #finding all x, y and z values in the output and joining each together
        # the x-line is always preceded by a nr and atomtype, so it works a bit different
        x_matches = re.findall(r'\d+.x\s\w{1,2}\s+-{0,1}\d\..+', summary)
        y_matches = re.findall(r'\s{5}y\s{6,7}-{0,1}\d\..+', summary)
        z_matches = re.findall(r'\s{5}z\s{6,7}-{0,1}\d\..+', summary)
        en_matches = re.findall(r'Energy\s+-?\d+\.\d+.*', summary)
        
        # Finding the number of points in the calculation, via the entire number of matches.
        irc_points = 0
        x_values, y_values, z_values, full_xyz, ens = [], [], [], [], []
        
        for i in range(len(en_matches)):
                ens.extend(string.split(en_matches[i])[1:])
                
        for i in range(len(ens)):
                ens[i] = float(ens[i])
        
        for i in range(len(z_matches)):
                irc_points += len(string.split(z_matches[i])[1:])
        irc_points = irc_points/len(atoms)      
        
        # making on list for x, y and z each.
        for i in range(len(atoms)):
                x_values.append(string.split(x_matches[i])[2:])
                y_values.append(string.split(y_matches[i])[1:])
                z_values.append(string.split(z_matches[i])[1:])
                for j in range(1, len(x_matches)/len(atoms)):
                        x_values[i].extend(string.split(x_matches[i+j*len(atoms)])[2:])
                        y_values[i].extend(string.split(y_matches[i+j*len(atoms)])[1:])
                        z_values[i].extend(string.split(z_matches[i+j*len(atoms)])[1:])
        
        #making the giant matrix
        for j in range(irc_points):  
                xyz_matrix = []
                for k in range(len(atoms)): 
                        xyz_matrix.append([float(x_values[k][j]), float(y_values[k][j]), float(z_values[k][j])])
                full_xyz.append(xyz_matrix)
        
        return atoms, full_xyz, ens

def average_distance(xyzmat):
    dist = [0]
    
    for i in range(1, len(xyzmat)):
        xyz1 = xyzmat[-1+i]
        xyz2 = xyzmat[i]
        diff = 0
        for j in range(len(xyz1)):
            diff += math.sqrt((xyz1[j][0] - xyz2[j][0])**2 + (xyz1[j][1] - xyz2[j][1])**2 + (xyz1[j][2] - xyz2[j][2])**2)
        dist.append(diff + dist[-1])
        
    return dist

def average_mw_distance(xyzmat, atoms):      
    dist = [0]
    for i in range(len(atoms)):
        atoms[i] = string.capitalize(atoms[i])
    
    am = {}
    am['C'] = 12.000
    am['H'] = 1.0078
    am['Pd'] = 105.9035
    am['Cl'] = 34.9689
    
    for i in range(1, len(xyzmat)):
        xyz1 = xyzmat[-1+i]
        xyz2 = xyzmat[i]      
        diff = 0
        for j in range(len(xyz1)):
            m = math.sqrt(am[atoms[j]])
            diff += math.sqrt((xyz1[j][0]*m - xyz2[j][0]*m)**2 + (xyz1[j][1]*m - xyz2[j][1]*m)**2 + (xyz1[j][2]*m - xyz2[j][2]*m)**2)
        dist.append(diff + dist[-1])

    return dist
   
def subs_geomvar(bondlength_vars, bondlengths, matrix_list):
    for i in range(len(bondlengths[0])):
        for j in range(len(bondlengths)):
            if float(bondlengths[j][i]) < 0 and re.search(r'\+'+bondlength_vars[j], matrix_list[i]): 
                matrix_list[i] = re.sub(r'\+'+bondlength_vars[j], str(bondlengths[j][i]), matrix_list[i])
            else:
                matrix_list[i] = re.sub(bondlength_vars[j], str(bondlengths[j][i]), matrix_list[i])
                
    return matrix_list
    

def shelljoin(*args):
        """ Sub backslashes with forwards for shell syntax"""
        path = ''
        for a in args:
                path = os.path.join(path, a)
        path = re.sub(r'\\', '/', path)
        return path

