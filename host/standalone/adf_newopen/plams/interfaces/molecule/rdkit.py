
__all__ = ['add_Hs', 'apply_reaction_smarts', 'apply_template',
           'gen_coords_rdmol', 'get_backbone_atoms', 'modify_atom',
           'to_rdmol', 'from_rdmol', 'from_sequence', 'from_smiles', 'from_smarts',
           'partition_protein', 'readpdb', 'writepdb', 'get_substructure', 'get_conformations',
           'yield_coords']

"""
@author: Lars Ridder
@description: A set of functions to manipulate molecules based on RDKit

This is a series of functions that apply RDKit functionality on PLAMS molecules
"""

import sys
import random
from warnings import warn
try:
    import dill as pickle
except ImportError:
    import pickle

try:
    from rdkit import Chem, Geometry
    from rdkit.Chem import AllChem
except ImportError:
    __all__ = []

from ...mol.bond import Bond
from ...mol.atom import Atom
from ...mol.molecule import Molecule
from ...core.functions import add_to_class


def from_rdmol(rdkit_mol, confid=-1, properties=True):
    """
    Translate an RDKit molecule into a PLAMS molecule type.
    RDKit properties will be unpickled if their name ends with '_pickled'.

    :parameter rdkit_mol: RDKit molecule
    :type rdkit_mol: rdkit.Chem.Mol
    :parameter int confid: conformer identifier from which to take coordinates
    :parameter bool properties: If all Chem.Mol, Chem.Atom and Chem.Bond properties should be converted from RDKit to PLAMS format.
    :return: a PLAMS molecule
    :rtype: |Molecule|
    """
    if isinstance(rdkit_mol, Molecule):
        return rdkit_mol
    # Create PLAMS molecule
    plams_mol = Molecule()
    total_charge = 0
    try:
        Chem.Kekulize(rdkit_mol)
    except Exception:
        pass
    conf = rdkit_mol.GetConformer(id=confid)

    # Add atoms and assign properties to the PLAMS atom if *properties* = True
    for rd_atom in rdkit_mol.GetAtoms():
        pos = conf.GetAtomPosition(rd_atom.GetIdx())
        ch = rd_atom.GetFormalCharge()
        pl_atom = Atom(
            rd_atom.GetAtomicNum(), coords=(pos.x, pos.y, pos.z), charge=ch)
        if properties and rd_atom.GetPDBResidueInfo():
            pl_atom.properties.pdb_info = get_PDBResidueInfo(rd_atom)
        plams_mol.add_atom(pl_atom)
        total_charge += ch

        # Check for R/S information
        stereo = str(rd_atom.GetChiralTag())
        if stereo == 'CHI_TETRAHEDRAL_CCW':
            pl_atom.properties.stereo = 'counter-clockwise'
        elif stereo == 'CHI_TETRAHEDRAL_CW':
            pl_atom.properties.stereo = 'clockwise'

    # Add bonds to the PLAMS molecule
    for bond in rdkit_mol.GetBonds():
        at1 = plams_mol.atoms[bond.GetBeginAtomIdx()]
        at2 = plams_mol.atoms[bond.GetEndAtomIdx()]
        plams_mol.add_bond(Bond(at1, at2, bond.GetBondTypeAsDouble()))

        # Check for cis/trans information
        stereo, bond_dir = str(bond.GetStereo()), str(bond.GetBondDir())
        if stereo == 'STEREOZ' or stereo == 'STEREOCIS':
            plams_mol.bonds[-1].properties.stereo = 'Z'
        elif stereo == 'STEREOE' or stereo == 'STEREOTRANS':
            plams_mol.bonds[-1].properties.stereo = 'E'
        elif bond_dir == 'ENDUPRIGHT':
            plams_mol.bonds[-1].properties.stereo = 'up'
        elif bond_dir == 'ENDDOWNRIGHT':
            plams_mol.bonds[-1].properties.stereo = 'down'

    # Set charge and assign properties to PLAMS molecule and bonds if *properties* = True
    plams_mol.charge = total_charge
    if properties:
        prop_from_rdmol(plams_mol, rdkit_mol)
        for rd_atom, plams_atom in zip(rdkit_mol.GetAtoms(), plams_mol):
            prop_from_rdmol(plams_atom, rd_atom)
        for rd_bond, plams_bond in zip(rdkit_mol.GetBonds(), plams_mol.bonds):
            prop_from_rdmol(plams_bond, rd_bond)
    return plams_mol


def to_rdmol(plams_mol, sanitize=True, properties=True, assignChirality=False):
    """
    Translate a PLAMS molecule into an RDKit molecule type.
    PLAMS |Molecule|, |Atom| or |Bond| properties are pickled if they are neither booleans, floats,
    integers, floats nor strings, the resulting property names are appended with '_pickled'.

    :parameter plams_mol: A PLAMS molecule
    :parameter bool sanitize: Kekulize, check valencies, set aromaticity, conjugation and hybridization
    :parameter bool properties: If all |Molecule|, |Atom| and |Bond| properties should be converted from PLAMS to RDKit format.
    :parameter bool assignChirality: Assign R/S and cis/trans information, insofar as this was not yet present in the PLAMS molecule.
    :type plams_mol: |Molecule|
    :return: an RDKit molecule
    :rtype: rdkit.Chem.Mol
    """
    if isinstance(plams_mol, Chem.Mol):
        return plams_mol
    # Create rdkit molecule
    e = Chem.EditableMol(Chem.Mol())

    # Add atoms and assign properties to the RDKit atom if *properties* = True
    for pl_atom in plams_mol.atoms:
        rd_atom = Chem.Atom(int(pl_atom.atnum))
        if 'charge' in pl_atom.properties:
            rd_atom.SetFormalCharge(pl_atom.properties.charge)
        if properties:
            if 'pdb_info' in pl_atom.properties:
                set_PDBresidueInfo(rd_atom, pl_atom.properties.pdb_info)
            for prop in pl_atom.properties:
                if prop not in ('charge', 'pdb_info', 'stereo'):
                    prop_to_rdmol(pl_atom, rd_atom, prop)

        # Check for R/S information
        if pl_atom.properties.stereo:
            stereo = pl_atom.properties.stereo.lower()
            if stereo == 'counter-clockwise':
                rd_atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
            elif stereo == 'clockwise':
                rd_atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
        e.AddAtom(rd_atom)

    # Mapping of PLAMS bond orders to RDKit bond types:
    def plams_to_rd_bonds(bo):
        if bo > 1.4 and bo < 1.6:
            return 12 # bond type for aromatic bond
        else:
            return int(bo)

    # Add bonds to the RDKit molecule
    for bond in plams_mol.bonds:
        a1 = plams_mol.atoms.index(bond.atom1)
        a2 = plams_mol.atoms.index(bond.atom2)
        e.AddBond(a1, a2, Chem.BondType(plams_to_rd_bonds(bond.order)))
    rdmol = e.GetMol()

    # Check for cis/trans information
    for pl_bond, rd_bond in zip(plams_mol.bonds, rdmol.GetBonds()):
        if pl_bond.properties.stereo:
            stereo = pl_bond.properties.stereo.lower()
            if stereo == 'e' or stereo == 'trans':
                rd_bond.SetStereo(Chem.rdchem.BondStereo.STEREOE)
            elif stereo == 'z' or stereo == 'cis':
                rd_bond.SetStereo(Chem.rdchem.BondStereo.STEREOZ)
            elif stereo == 'up':
                rd_bond.SetBondDir(Chem.rdchem.BondDir.ENDUPRIGHT)
            elif stereo == 'down':
                rd_bond.SetBondDir(Chem.rdchem.BondDir.ENDDOWNRIGHT)

    # Assign properties to RDKit molecule and bonds if *properties* = True
    if properties:
        for prop in plams_mol.properties:
            prop_to_rdmol(plams_mol, rdmol, prop)
        for pl_bond, rd_bond in zip(plams_mol.bonds, rdmol.GetBonds()):
            for prop in pl_bond.properties:
                if prop != 'stereo':
                    prop_to_rdmol(pl_bond, rd_bond, prop)

    if sanitize:
        Chem.SanitizeMol(rdmol)
    conf = Chem.Conformer()
    for i, atom in enumerate(plams_mol.atoms):
        xyz = Geometry.Point3D(atom._getx(), atom._gety(), atom._getz())
        conf.SetAtomPosition(i, xyz)
    rdmol.AddConformer(conf)
    # REB: Assign all stereochemistry, if it wasn't already there
    if assignChirality :
        Chem.rdmolops.AssignAtomChiralTagsFromStructure(rdmol,confId=conf.GetId(),replaceExistingTags=False)
        try :
            Chem.AssignStereochemistryFrom3D(rdmol,confId=conf.GetId(),replaceExistingTags=False)
        except AttributeError :
            pass
    return rdmol


pdb_residue_info_items = [
    'AltLoc', 'ChainId', 'InsertionCode', 'IsHeteroAtom', 'Name', 'Occupancy',
    'ResidueName', 'ResidueNumber', 'SecondaryStructure', 'SegmentNumber',
    'SerialNumber', 'TempFactor']
# 'MonomerType' was excluded because it is an rdkit type that cannot easilty be serialized


def get_PDBResidueInfo(rdkit_atom):
    pdb_info = {}
    for item in pdb_residue_info_items:
        get_function = 'Get' + item
        pdb_info[item] = rdkit_atom.GetPDBResidueInfo().__getattribute__(get_function)()
    return pdb_info


def set_PDBresidueInfo(rdkit_atom, pdb_info):
    atom_pdb_residue_info = Chem.AtomPDBResidueInfo()
    for item, value in pdb_info.items():
        set_function = 'Set' + item
        atom_pdb_residue_info.__getattribute__(set_function)(value)
    rdkit_atom.SetMonomerInfo(atom_pdb_residue_info)


def prop_to_rdmol(pl_obj, rd_obj, prop):
    """
    Convert a single PLAMS property into an RDKit property.

    :paramter pl_obj: A PLAMS object.
    :type pl_obj: |Molecule|, |Atom| or |Bond|.
    :parameter rd_obj: An RDKit object.
    :type rd_obj: rdkit.Chem.Mol, rdkit.Chem.Atom or rdkit.Chem.Bond
    :parameter str prop: The |Settings| key of the PLAMS property.
    """
    obj = type(pl_obj.properties.get(prop))
    obj_dict = {bool: rd_obj.SetBoolProp,
                float: rd_obj.SetDoubleProp,
                int: rd_obj.SetIntProp,
                str: rd_obj.SetProp}
    if obj_dict.get(obj):
        obj_dict[obj](prop, pl_obj.properties.get(prop))
    else:
        name = prop + '_pickled'
        try:
            rd_obj.SetProp(name, pickle.dumps(pl_obj.properties.get(prop), 0).decode())
        except (Exception, pickle.PicklingError):
            pass


def prop_from_rdmol(pl_obj, rd_obj):
    """
    Convert one or more RDKit properties into PLAMS properties.

    :paramter pl_obj: A PLAMS object.
    :type pl_obj: |Molecule|, |Atom| or |Bond|.
    :parameter rd_obj: An RDKit object.
    :type rd_obj: rdkit.Chem.Mol, rdkit.Chem.Atom or rdkit.Chem.Bond
    """
    prop_dict = rd_obj.GetPropsAsDict()
    for propname in prop_dict.keys():
        if propname == '__computedProps':
            continue
        if '_pickled' not in propname:
            pl_obj.properties[propname] = prop_dict[propname]
        else:
            prop = prop_dict[propname]
            propname = propname.rsplit('_pickled', 1)[0]
            pl_obj.properties[propname] = pickle.loads(prop.encode())


def from_smiles(smiles, nconfs=1, name=None, forcefield=None, rms=0.1):
    """
    Generates PLAMS molecule(s) from a smiles strings.

    :parameter str smiles: A smiles string
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry optimization
        and ranking of comformations. The default value None results in skipping of the
        geometry optimization step.
    :parameter float rms: Root Mean Square deviation threshold for
        removing similar/equivalent conformations
    :return: A molecule with hydrogens and 3D coordinates or a list of molecules if nconfs > 1
    :rtype: |Molecule| or list of PLAMS Molecules
    """
    smiles = str(smiles.split()[0])
    smiles = Chem.CanonSmiles(smiles)
    rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    rdkit_mol.SetProp('smiles', smiles)
    return get_conformations(rdkit_mol, nconfs, name, forcefield, rms)


def from_smarts(smarts, nconfs=1, name=None, forcefield=None, rms=0.1):
    """
    Generates PLAMS molecule(s) from a smarts strings.
    This allows for example to define hydrogens explicitly.
    However it is less suitable for aromatic molecules (use from_smiles in that case).

    :parameter str smarts: A smarts string
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry
        optimization and ranking of comformations. The default value None results
        in skipping of the geometry optimization step.
    :parameter float rms: Root Mean Square deviation threshold for removing
        similar/equivalent conformations.
    :return: A molecule with hydrogens and 3D coordinates or a list of molecules if nconfs > 1
    :rtype: |Molecule| or list of PLAMS Molecules
    """
    smiles = str(smarts.split()[0])
    mol = Chem.MolFromSmarts(smiles)
    Chem.SanitizeMol(mol)
    molecule = Chem.AddHs(mol)
    molecule.SetProp('smiles', smiles)
    return get_conformations(molecule, nconfs, name, forcefield, rms)


def get_conformations(mol, nconfs=1, name=None, forcefield=None, rms=-1, enforceChirality=False):
    """
    Generates 3D conformation(s) for an rdkit_mol or a PLAMS Molecule

    :parameter mol: RDKit or PLAMS Molecule
    :type mol: rdkit.Chem.Mol or |Molecule|
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry
        optimization and ranking of comformations. The default value None results
        in skipping of the geometry optimization step
    :parameter float rms: Root Mean Square deviation threshold for removing
        similar/equivalent conformations.
    :parameter bool enforceChirality: Enforce the correct chirality if chiral centers are present
    :return: A molecule with hydrogens and 3D coordinates or a list of molecules if nconfs > 1
    :rtype: |Molecule| or list of PLAMS Molecules
    """

    if isinstance(mol, Molecule):
        rdkit_mol = to_rdmol(mol,assignChirality=enforceChirality)
    else:
        rdkit_mol = mol

    def MMFFenergy(cid):
        ff = AllChem.MMFFGetMoleculeForceField(
            rdkit_mol, AllChem.MMFFGetMoleculeProperties(rdkit_mol), confId=cid)
        try:
            energy = ff.CalcEnergy()
        except:
            msg = "MMFF energy calculation failed for molecule: " + Chem.MolToSmiles(rdkit_mol) + \
                  "\nNo geometry optimization was performed."
            warn(msg)
            energy = 1e9
        return energy

    def UFFenergy(cid):
        ff = AllChem.UFFGetMoleculeForceField(rdkit_mol, confId=cid)
        try:
            energy = ff.CalcEnergy()
        except:
            msg = "MMFF energy calculation failed for molecule: " + Chem.MolToSmiles(rdkit_mol) + \
                  "\nNo geometry optimization was performed."
            warn(msg)
            energy = 1e9
        return energy

    if name:
        rdkit_mol.SetProp('name', name)

    try:
        cids = list(AllChem.EmbedMultipleConfs(rdkit_mol,nconfs,pruneRmsThresh=rms,randomSeed=1,enforceChirality=enforceChirality))
    except:
         # ``useRandomCoords = True`` prevents (poorly documented) crash for large systems
        cids = list(AllChem.EmbedMultipleConfs(rdkit_mol,nconfs,pruneRmsThresh=rms,randomSeed=1,useRandomCoords=True,enforceChirality=enforceChirality))

    if forcefield:
        # Select the forcefield (UFF or MMFF)
        optimize_molecule, energy = {
            'uff': [AllChem.UFFOptimizeMolecule, UFFenergy],
            'mmff': [AllChem.MMFFOptimizeMolecule, MMFFenergy],
        }[forcefield]

        # Optimize and sort conformations
        for cid in cids:
            optimize_molecule(rdkit_mol, confId=cid)
        cids.sort(key=energy)

        # Remove duplicate conformations based on RMS
        if rms > 0:
            keep = [cids[0]]
            for cid in cids[1:]:
                for idx in keep:
                    try:
                        r = AllChem.AlignMol(rdkit_mol, rdkit_mol, cid, idx)
                    except:
                        r = rms + 1
                        message = "Alignment failed in multiple conformation generation: "
                        message += Chem.MolToSmiles(rdkit_mol)
                        message += "\nAssuming different conformations."
                        warn(message)
                    if r < rms:
                        break
                else:
                    keep.append(cid)
            cids = keep

    if nconfs == 1:
        return from_rdmol(rdkit_mol)
    else:
        return [from_rdmol(rdkit_mol, cid) for cid in cids]


def from_sequence(sequence, nconfs=1, name=None, forcefield=None, rms=0.1):
    """
    Generates PLAMS molecule from a peptide sequence.
    Includes explicit hydrogens and 3D coordinates.

    :parameter str sequence: A peptide sequence, e.g. 'HAG'
    :parameter int nconfs: Number of conformers to be generated
    :parameter str name: A name for the molecule
    :parameter str forcefield: Choose 'uff' or 'mmff' forcefield for geometry
        optimization and ranking of comformations. The default value None results
        in skipping of the geometry optimization step.
    :parameter float rms: Root Mean Square deviation threshold for removing
        similar/equivalent conformations.
    :return: A peptide molecule with hydrogens and 3D coordinates
        or a list of molecules if nconfs > 1
    :rtype: |Molecule| or list of PLAMS Molecules
    """
    rdkit_mol = Chem.AddHs(Chem.MolFromSequence(sequence))
    rdkit_mol.SetProp('sequence', sequence)
    return get_conformations(rdkit_mol, nconfs, name, forcefield, rms)


def calc_rmsd(mol1, mol2):
    """
    Superimpose two molecules and calculate the root-mean-squared deviations of
    the atomic positions.
    The molecules should be identical, but the ordering of the atoms may differ.

    :param mol1: Molecule 1
    :param mol2: Molecule 2
    :return: The rmsd after superposition
    :rtype: float
    """
    rdkit_mol1 = to_rdmol(mol1)
    rdkit_mol2 = to_rdmol(mol2)
    try:
        return AllChem.GetBestRMS(rdkit_mol1, rdkit_mol2)
    except:
        return -999


def modify_atom(mol, idx, element):
    """
    Change atom "idx" in molecule "mol" to "element" and add or remove hydrogens accordingly

    :parameter mol: molecule to be modified
    :type mol: |Molecule| or rdkit.Chem.Mol
    :parameter int idx: index of the atom to be modified
    :parameter str element:
    :return: Molecule with new element and possibly added or removed hydrogens
    :rtype: |Molecule|
    """
    rdmol = to_rdmol(mol)
    if rdmol.GetAtomWithIdx(idx).GetSymbol() == element:
        return mol
    else:
        e = Chem.EditableMol(rdmol)
        for neighbor in reversed(rdmol.GetAtomWithIdx(idx - 1).GetNeighbors()):
            if neighbor.GetSymbol() == 'H':
                e.RemoveAtom(neighbor.GetIdx())
        e.ReplaceAtom(idx - 1, Chem.Atom(element))
        newmol = e.GetMol()
        Chem.SanitizeMol(newmol)
        newmol = Chem.AddHs(newmol, addCoords=True)
        return from_rdmol(newmol)


def apply_template(mol, template):
    """
    Modifies bond orders in PLAMS molecule according template smiles structure.

    :parameter mol: molecule to be modified
    :type mol: |Molecule| or rdkit.Chem.Mol
    :parameter str template: smiles string defining the correct chemical structure
    :return: Molecule with correct chemical structure and provided 3D coordinates
    :rtype: |Molecule|
    """
    rdmol = to_rdmol(mol, sanitize=False)
    template_mol = Chem.AddHs(Chem.MolFromSmiles(template))
    newmol = Chem.AllChem.AssignBondOrdersFromTemplate(template_mol, rdmol)
    return from_rdmol(newmol)


def apply_reaction_smarts(
        mol, reaction_smarts, complete=False, forcefield=None, return_rdmol=False):
    """
    Applies reaction smirks and returns product.
    If returned as a PLAMS molecule, thismolecule.properties.orig_atoms
    is a list of indices of atoms that have not been changed
    (which can for example be used partially optimize new atoms only with the freeze keyword)

    :parameter mol: molecule to be modified
    :type mol: |Molecule| or rdkit.Chem.Mol
    :parameter str reactions_smarts: Reactions smarts to be applied to molecule
    :parameter complete: Apply reaction until no further changes occur or given
        fraction of reaction centers have been modified
    :type complete: bool or float (value between 0 and 1)
    :parameter forcefield: Specify 'uff' or 'mmff' to apply forcefield based
        geometry optimization of product structures.
    :type forcefield: str
    :param bool return_rdmol: return a RDKit molecule if true, otherwise a PLAMS molecule
    :return: (product molecule, list of unchanged atoms)
    :rtype: (|Molecule|, list of int)
    """
    def react(reactant, reaction):
        """ Apply reaction to reactant and return products """
        ps = reaction.RunReactants([reactant])
        # if reaction doesn't apply, return the reactant
        if len(ps) == 0:
            return [(reactant, range(reactant.GetNumAtoms()))]
        full = len(ps)
        while complete:  # when complete is True
            # apply reaction until no further changes
            r = random.randint(0, len(ps) - 1)
            reactant = ps[r][0]
            ps = reaction.RunReactants([reactant])
            if len(ps) == 0 or len(ps) / full < (1 - complete):
                ps = [[reactant]]
                break
        # add hydrogens and generate coordinates for new atoms
        products = []
        for p in ps[0]:
            Chem.SanitizeMol(p)
            q = Chem.AddHs(p)
            Chem.SanitizeMol(q)
            u = gen_coords_rdmol(q)  # These are the atoms that have not changed
            products.append((q, u))
        return products

    mol = to_rdmol(mol)
    reaction = AllChem.ReactionFromSmarts(reaction_smarts)
    # RDKit removes fragments that are disconnected from the reaction center
    # In order to keep these, the molecule is first split in separate fragments
    # and the results, including non-reacting parts, are re-combined afterwards
    frags = (Chem.GetMolFrags(mol, asMols=True))
    product = Chem.Mol()
    unchanged = []  # List of atoms that have not changed
    for frag in frags:
        for p, u in react(frag, reaction):
            unchanged += [product.GetNumAtoms() + i for i in u]
            product = Chem.CombineMols(product, p)
    if forcefield:
        optimize_coordinates(product, forcefield, fixed=unchanged)
    # The molecule is returned together with a list of atom indices of the atoms
    # that are identical to those
    # in the reactants. This list can be used in subsequent partial optimization of the molecule
    if not return_rdmol:
        product = from_rdmol(product)
        product.properties.orig_atoms = [a + 1 for a in unchanged]
    return product


def gen_coords(plamsmol):
    """ Calculate 3D positions only for atoms without coordinates """
    rdmol = to_rdmol(plamsmol)
    unchanged = gen_coords_rdmol(rdmol)
    conf = rdmol.GetConformer()
    for a in range(len(plamsmol.atoms)):
        pos = conf.GetAtomPosition(a)
        atom = plamsmol.atoms[a]
        atom._setx(pos.x)
        atom._sety(pos.y)
        atom._setz(pos.z)
    return [a + 1 for a in unchanged]


def gen_coords_rdmol(rdmol):
    ref = rdmol.__copy__()
    conf = rdmol.GetConformer()
    coordDict = {}
    unchanged = []
    maps = []
    # Put known coordinates in coordDict
    for i in range(rdmol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        if (-0.0001 < pos.x < 0.0001) and (-0.0001 < pos.y < 0.0001) and \
           (-0.0001 < pos.z < 0.0001):
            continue  # atom without coordinates
        coordDict[i] = pos
        unchanged.append(i)
        maps.append((i, i))
    # compute coordinates for new atoms, keeping known coordinates
    rms = 1
    rs = 1
    # repeat embedding and alignment until the rms of mapped atoms is sufficiently small
    if rdmol.GetNumAtoms() > len(maps):
        while rms > 0.1:
            AllChem.EmbedMolecule(rdmol, coordMap=coordDict, randomSeed=rs,
                                  useBasicKnowledge=True)
            # align new molecule to original coordinates
            rms = AllChem.AlignMol(rdmol, ref, atomMap=maps)
            rs += 1
    return unchanged


def optimize_coordinates(rdkit_mol, forcefield, fixed=[]):
    def MMFFminimize():
        ff = AllChem.MMFFGetMoleculeForceField(
            rdkit_mol, AllChem.MMFFGetMoleculeProperties(rdkit_mol))
        for f in fixed:
            ff.AddFixedPoint(f)
        try:
            ff.Minimize()
        except:
            warn("MMFF geometry optimization failed for molecule: " + Chem.MolToSmiles(rdkit_mol))

    def UFFminimize():
        ff = AllChem.UFFGetMoleculeForceField(rdkit_mol, ignoreInterfragInteractions=True)
        for f in fixed:
            ff.AddFixedPoint(f)
        try:
            ff.Minimize()
        except:
            warn("UFF geometry optimization failed for molecule: " + Chem.MolToSmiles(rdkit_mol))
    optimize_molecule = {
        'uff': UFFminimize,
        'mmff': MMFFminimize}[forcefield]
    Chem.SanitizeMol(rdkit_mol)
    optimize_molecule()
    return


def write_molblock(plams_mol, file=sys.stdout):
    file.write(Chem.MolToMolBlock(to_rdmol(plams_mol)))


def readpdb(pdb_file, removeHs=False, proximityBonding=False, return_rdmol=False):
    """
    Generate a molecule from a PDB file

    :param pdb_file: The PDB file to read
    :type pdb_file: path- or file-like
    :param bool removeHs: Hydrogens are removed if True
    :param bool proximityBonding: Enables automatic proximity bonding
    :param bool return_rdmol: return a RDKit molecule if true, otherwise a PLAMS molecule
    :return: The molecule
    :rtype: |Molecule| or rdkit.Chem.Mol
    """
    try:
        pdb_file = open(pdb_file, 'r')
    except TypeError:
        pass  # pdb_file is a file-like object... hopefully

    pdb_mol = Chem.MolFromPDBBlock(pdb_file.read(), removeHs=removeHs, proximityBonding=proximityBonding)
    return pdb_mol if return_rdmol else from_rdmol(pdb_mol)


def writepdb(mol, pdb_file=sys.stdout):
    """
    Write a PDB file from a molecule

    :parameter mol: molecule to be exported to PDB
    :type mol: |Molecule| or rdkit.Chem.Mol
    :param pdb_file: The PDB file to write to, or a filename
    :type pdb_file: path- or file-like
    """
    try:
        pdb_file = open(pdb_file, 'w')
    except TypeError:
        pass  # pdb_file is a file-like object... hopefully

    mol = to_rdmol(mol)
    pdb_file.write(Chem.MolToPDBBlock(mol))


def add_Hs(mol, forcefield=None, return_rdmol=False):
    """
    Add hydrogens to protein molecules read from PDB.
    Makes sure that the hydrogens get the correct PDBResidue info.

    :param mol: Molecule to be protonated
    :type mol: |Molecule| or rdkit.Chem.Mol
    :param str forcefield: Specify 'uff' or 'mmff' to apply forcefield based
        geometry optimization on new atoms.
    :param bool return_rdmol: return a RDKit molecule if true, otherwise a PLAMS molecule
    :return: A molecule with explicit hydrogens added
    :rtype: |Molecule| or rdkit.Chem.Mol
    """
    mol = to_rdmol(mol)
    retmol = Chem.AddHs(mol)
    for atom in retmol.GetAtoms():
        if atom.GetPDBResidueInfo() is None and atom.GetSymbol() == 'H':
            bond = atom.GetBonds()[0]
            if bond.GetBeginAtom().GetIdx() == atom.GetIdx:
                connected_atom = bond.GetEndAtom()
            else:
                connected_atom = bond.GetBeginAtom()
            try:
                ResInfo = connected_atom.GetPDBResidueInfo()
                if ResInfo is None:
                    continue  # Segmentation faults are raised if ResInfo is None
                atom.SetMonomerInfo(ResInfo)
                atomname = 'H' + atom.GetPDBResidueInfo().GetName()[1:]
                atom.GetPDBResidueInfo().SetName(atomname)
            except:
                pass
    unchanged = gen_coords_rdmol(retmol)
    if forcefield:
        optimize_coordinates(retmol, forcefield, fixed=unchanged)
    return retmol if return_rdmol else from_rdmol(retmol)


def add_fragment(rwmol, frag, rwmol_atom_idx=None, frag_atom_idx=None,
                 bond_order=None):
    molconf = rwmol.GetConformer()
    fragconf = frag.GetConformer()
    new_indices = []
    for a in frag.GetAtoms():
        new_index = rwmol.AddAtom(a)
        new_indices.append(new_index)
        molconf.SetAtomPosition(new_index, fragconf.GetAtomPosition(a.GetIdx()))
    for b in frag.GetBonds():
        ba = b.GetBeginAtomIdx()
        ea = b.GetEndAtomIdx()
        rwmol.AddBond(new_indices[ba], new_indices[ea], b.GetBondType())
    if bond_order:
        rwmol.AddBond(rwmol_atom_idx, new_indices[frag_atom_idx],
                      Chem.BondType.values[bond_order])
        rwmol.GetAtomWithIdx(new_indices[frag_atom_idx]).SetNumRadicalElectrons(0)


def get_fragment(mol, indices, incl_expl_Hs=True, neutralize=True):
    molconf = mol.GetConformer()
    fragment = Chem.RWMol(Chem.Mol())
    fragconf = Chem.Conformer()
    # Put atoms in fragment
    for i in indices:
        atom = mol.GetAtomWithIdx(i)
        new_index = fragment.AddAtom(atom)
        pos = molconf.GetAtomPosition(i)
        fragconf.SetAtomPosition(new_index, pos)
    # Put bonds in fragment
    for b in mol.GetBonds():
        ba = b.GetBeginAtomIdx()
        ea = b.GetEndAtomIdx()
        if ba in indices and ea in indices:
            fragment.AddBond(indices.index(ba), indices.index(ea),
                             b.GetBondType())
            continue
        if not incl_expl_Hs:
            continue
        if ba in indices and mol.GetAtomWithIdx(ea).GetSymbol() == 'H':
            hi = fragment.AddAtom(mol.GetAtomWithIdx(ea))
            fragconf.SetAtomPosition(hi, molconf.GetAtomPosition(ea))
            fragment.AddBond(indices.index(ba), hi, Chem.BondType.SINGLE)
            continue
        if ea in indices and mol.GetAtomWithIdx(ba).GetSymbol() == 'H':
            hi = fragment.AddAtom(mol.GetAtomWithIdx(ba))
            fragconf.SetAtomPosition(hi, molconf.GetAtomPosition(ba))
            fragment.AddBond(indices.index(ea), hi, Chem.BondType.SINGLE)
    ret_frag = fragment.GetMol()
    Chem.SanitizeMol(ret_frag)
    if neutralize:
        for atom in ret_frag.GetAtoms():
            nrad = atom.GetNumRadicalElectrons()
            if nrad > 0:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + nrad)
                atom.SetNumRadicalElectrons(0)
    Chem.SanitizeMol(ret_frag)
    ret_frag.AddConformer(fragconf)
    return ret_frag


def partition_protein(mol, residue_bonds=None, split_heteroatoms=True, return_rdmol=False):
    """
    Splits a protein molecule into capped amino acid fragments and caps.

    :param mol: A protein molecule
    :type mol: |Molecule| or rdkit.Chem.Mol
    :param tuple residue_bonds: a tuple of pairs of residue number indicating which
        peptide bonds to split. If none, split all peptide bonds.
    :param bool split_heteroatoms: if True, all bonds between a heteroatom and
        a non-heteroatom across residues are removed
    :return: list of fragments, list of caps
    """
    mol = to_rdmol(mol)
    caps = []
    em = Chem.RWMol(mol)
    if split_heteroatoms:
        for bond in mol.GetBonds():
            resinfa = bond.GetBeginAtom().GetPDBResidueInfo()
            resinfb = bond.GetEndAtom().GetPDBResidueInfo()
            if resinfa.GetIsHeteroAtom() is not resinfb.GetIsHeteroAtom():
                if resinfa.GetResidueNumber() != resinfb.GetResidueNumber():
                    em.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    # Split peptide bonds
    pept_bond = Chem.MolFromSmarts('[C;X4;H1,H2][CX3](=O)[NX3][C;X4;H1,H2][CX3](=O)')
    for match in mol.GetSubstructMatches(pept_bond):
        if residue_bonds:
            resa = mol.GetAtomWithIdx(match[1]).GetPDBResidueInfo().GetResidueNumber()
            resb = mol.GetAtomWithIdx(match[3]).GetPDBResidueInfo().GetResidueNumber()
            if (resa, resb) not in residue_bonds and (resb, resa) not in residue_bonds:
                continue
        cap = get_fragment(mol, match[0:5])
        cap = add_Hs(cap, return_rdmol=True)
        caps.append(cap if return_rdmol else from_rdmol(cap))
        cap_o_ind = cap.GetSubstructMatch(Chem.MolFromSmarts('[C;X4][CX3]=O'))
        cap_o = get_fragment(cap, cap_o_ind, neutralize=False)
        cap_n_ind = cap.GetSubstructMatch(Chem.MolFromSmarts('O=[CX3][NX3][C;X4]'))[2:]
        cap_n = get_fragment(cap, cap_n_ind, neutralize=False)
        em.RemoveBond(match[1], match[3])
        add_fragment(em, cap_o, match[3], 1, 1)
        add_fragment(em, cap_n, match[1], 0, 1)
    # Split disulfide bonds
    ss_bond = Chem.MolFromSmarts('[C;X4;H1,H2]SS[C;X4;H1,H2]')
    for match in mol.GetSubstructMatches(ss_bond):
        cap = get_fragment(mol, match[0:5])
        cap = add_Hs(cap, return_rdmol=True)
        caps.append(cap if return_rdmol else from_rdmol(cap))
        cap_s_ind = cap.GetSubstructMatch(Chem.MolFromSmarts('[C;X4]SS[C;X4]'))
        cap_s1 = get_fragment(cap, cap_s_ind[0:2], neutralize=False)
        cap_s2 = get_fragment(cap, cap_s_ind[2:4], neutralize=False)
        em.RemoveBond(match[1], match[2])
        add_fragment(em, cap_s1, match[2], 1, 1)
        add_fragment(em, cap_s2, match[1], 0, 1)
    frags = Chem.GetMolFrags(em.GetMol(), asMols=True, sanitizeFrags=False)
    if not return_rdmol:
        frags = [from_rdmol(frag) for frag in frags]
    return frags, caps


def charge_AAs(mol, return_rdmol=False):
    ionizations = {
        'ARG_NH2': 1,
        'LYS_NZ': 1,
        'GLU_OE2': -1,
        'ASP_OD2': -1}
    mol = to_rdmol(mol)
    for atom in mol.GetAtoms():
        resinfo = atom.GetPDBResidueInfo()
        res_atom = resinfo.GetResidueName() + '_' + resinfo.GetName().strip()
        try:
            atom.SetFormalCharge(ionizations[res_atom])
            Chem.SanitizeMol(mol)
        except KeyError:
            pass
        Chem.SanitizeMol(mol)
    return mol if return_rdmol else from_rdmol(mol)


def get_backbone_atoms(mol):
    """
    Return a list of atom indices corresponding to the backbone atoms in a peptide molecule.
    This function assumes PDB information in properties.pdb_info of each atom, which is the case
    if the molecule is generated with the "readpdb" or "from_sequence" functions.

    :parameter mol: a peptide molecule
    :type mol: |Molecule| or rdkit.Chem.Mol
    :return: a list of atom indices
    :rtype: list
    """
    mol = from_rdmol(mol)
    backbone = ['N', 'CA', 'C', 'O']
    return [a for a in range(1, len(mol) + 1)
            if str(mol[a].properties.pdb_info.Name).strip() in backbone]


def get_substructure(mol, func_list):
    """
    Search for functional groups within a molecule based on a list of reference functional groups.
    SMILES strings, PLAMS and/or RDKit molecules can be used interchangeably in "func_list".

    Example:

    .. code:: python

        >>> mol = from_smiles('OCCO')  # Ethylene glycol
        >>> func_list = ['[H]O', 'C[N+]', 'O=PO']
        >>> get_substructure(mol, func_list)

        {'[H]O': [(<scm.plams.mol.atom.Atom at 0x125183518>,
                   <scm.plams.mol.atom.Atom at 0x1251836a0>),
                  (<scm.plams.mol.atom.Atom at 0x125183550>,
                   <scm.plams.mol.atom.Atom at 0x125183240>)]}

    :parameter mol: A PLAMS molecule.
    :type mol: |Molecule|
    :parameter list func_list: A list of functional groups.
        Functional groups can be represented by SMILES strings, PLAMS and/or RDKit molecules.
    :return: A dictionary with functional groups from "func_list" as keys and a list of n-tuples
        with matching PLAMS |Atom| as values.
    """
    def _to_rdmol(functional_group):
        """ Turn a SMILES strings, RDKit or PLAMS molecules into an RDKit molecule. """
        if isinstance(functional_group, str):
            # RDKit tends to remove explicit hydrogens if SANITIZE_ADJUSTHS is enabled
            sanitize = Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_ADJUSTHS
            ret = Chem.MolFromSmiles(functional_group, sanitize=False)
            Chem.rdmolops.SanitizeMol(ret, sanitizeOps=sanitize)
            return ret
        elif isinstance(functional_group, Molecule):
            return to_rdmol(functional_group)
        elif isinstance(functional_group, Chem.Mol):
            return functional_group
        raise TypeError('get_substructure: ' + str(type(functional_group)) + ' is not a supported \
                        object type')

    def _get_match(mol, rdmol, functional_group):
        """ Perform a substructure match on "mol".
        If a match is found, return a list of n-tuples consisting PLAMS |Atom|.
        Otherwise return False. """
        matches = rdmol.GetSubstructMatches(functional_group)
        if matches:
            return [tuple(mol[j+1] for j in idx_tup) for idx_tup in matches]
        return False

    rdmol = to_rdmol(mol)
    rdmol_func_list = [_to_rdmol(i) for i in func_list]
    gen = (_get_match(mol, rdmol, i) for i in rdmol_func_list)
    return {key: value for key, value in zip(func_list, gen) if value}

  
def yield_coords(rdmol, id=-1):
    """Take an rdkit molecule and yield its coordinates as 3-tuples.

    .. code-block:: python

        >>> from scm.plams import yield_coords
        >>> from rdkit import Chem

        >>> rdmol = Chem.Mol(...)  # e.g. Methane
        >>> for xyz in yield_coords(rdmol):
        ...     print(xyz)
        (-0.0, -0.0, -0.0)
        (0.6405, 0.6405, -0.6405)
        (0.6405, -0.6405, 0.6405)
        (-0.6405, 0.6405, 0.6405)
        (-0.6405, -0.6405, -0.6405)


    The iterator produced by this function can, for example, be passed to
    :meth:`Molecule.from_array()<scm.plams.mol.molecule.Molecule.from_array>`
    the update the coordinates of a PLAMS Molecule in-place.

    .. code-block:: python

        >>> from scm.plams import Molecule

        >>> mol = Molecule(...)

        >>> xyz_iterator = yield_coords(rdmol)
        >>> mol.from_array(xyz_iterator)


    :parameter rdmol: An RDKit mol.
    :type rdmol: rdkit.Chem.Mol
    :parameter int id: The ID of the desired conformer.
    :return: An iterator yielding 3-tuples with *rdmol*'s Cartesian coordinates.
    :rtype: iterator
    """
    conf = rdmol.GetConformer(id=id)
    for atom in rdmol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        yield (pos.x, pos.y, pos.z)

        
@add_to_class(Molecule)
def assign_chirality(self):
    """
    Assigns stereo-info to PLAMS molecule by invoking RDKIT
    """
    rd_mol = to_rdmol(self, assignChirality=True)
    pl_mol = from_rdmol(rd_mol)

    # Add R/S info to self
    for iat,pl_atom in enumerate(pl_mol.atoms):
        # Check for R/S information
        if pl_atom.properties.stereo:
            self.atoms[iat].properties.stereo = pl_atom.properties.stereo

    # Add cis/trans information to self
    for ibond,pl_bond in enumerate(pl_mol.bonds):
        if pl_bond.properties.stereo:
            self.bonds[ibond] = pl_bond.properties.stereo

            
@add_to_class(Molecule)
def get_chirality(self):
    """
    Returns the chirality of the atoms
    """
    rd_mol = to_rdmol(self, assignChirality=True)
    return Chem.FindMolChiralCenters(rd_mol,force=True,includeUnassigned=True)
