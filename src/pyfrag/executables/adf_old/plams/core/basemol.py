from __future__ import unicode_literals

import copy
import heapq
import math
import numpy
import os
import types
import uuid

from .common import log
from .errors import MoleculeError, PTError, FileError
from .settings import Settings
from ..tools.pdbtools import PDBHandler, PDBRecord
from ..tools.utils import Units, PT

__all__ = ['Atom', 'Bond', 'Molecule']



#===================================================================================================
#===================================================================================================
#===================================================================================================



class Atom(object):
    """A class representing a single atom in three dimensional space.

    An instance of this class has the following attributes:
        *   ``atnum`` -- atomic number (zero for "dummy atoms")
        *   ``coords`` -- tuple of length 3 storing spatial coordinates
        *   ``bonds`` -- list of bonds (see |Bond|) this atom is a part of
        *   ``mol`` -- a |Molecule| this atom belongs to
        *   ``properties`` -- a |Settings| instance storing all other information about this atom (initially it is populated with *\*\*other* keyword arguments passed to the constructor)

    All the above attributes can be accessed either directly or using one of the following properties:
        *   ``x``, ``y``, ``z`` -- allow to read or modify each coordinate separately
        *   ``symbol`` -- allows to read or write atomic symbol directly. Atomic symbol is not stored as an attribute, instead of that atomic number (``atnum``) indicates the type of atom. In fact, ``symbol`` this is just a wrapper around ``atnum`` that uses |PeriodicTable| as a translator::

                >>> a = Atom(atnum=8)
                >>> print a.symbol
                O
                >>> a.symbol = 'Ca'
                >>> print a.atnum
                20

        *   ``mass`` -- atomic mass, obtained from |PeriodicTable|, read only
        *   ``radius`` -- atomic radius, obtained from |PeriodicTable|, read only
        *   ``connectors`` -- number of connectors, obtained from |PeriodicTable|, read only

    .. note::

        When creating a new atom, its type can be chosen either by setting an atomic number or a symbol (``atnum`` and ``symbol`` constructor arguments). Symbol takes precedence -- if it is supplied, ``atnum`` argument is ignored.


    Values stored in ``coords`` tuple do not necessarily have to be numeric, you can also store any string there. This might come handy for programs that allow parametrization of coordinates in the input file (to enforce some geometry constraints for example)::

            >>> a = Atom(symbol='C', coords=(1,2,3))
            >>> print a
                     C       1.00000       2.00000       3.00000
            >>> a.y = 'param1'
            >>> print a
                     C       1.00000        param1       3.00000

    However, non-numerical coordinates cannot be used together with some methods (for example :meth:`distance_to` or :meth:`translate`). Trying to do this will raise an exception.

    Internally, atomic coordinates are always expressed in angstroms. Most of methods that read or modify atomic coordinates accept keyword argument ``unit`` allowing to choose unit in which results and/or arguments are expressed (see |Units| for details). Throughout the entire code angstrom is the default length unit. If you don't specify ``unit`` parameter in any place of your script, all automatic unit handling described above boils down to occasional multiplication/division by 1.0.
    """
    def __init__(self, atnum=0, symbol=None, coords=None, unit='angstrom', bonds=None, mol=None, **other):
        if symbol is not None:
            self.symbol = symbol
        else:
            self.atnum = atnum
        self.mol = mol
        self.bonds = bonds or []
        self.properties = Settings(other)

        if coords is None:
            self.coords = (0.0, 0.0, 0.0)
        elif len(coords) == 3:
            tmp = []
            for i in coords:
                try:
                    i = Units.convert(float(i), unit, 'angstrom')
                except ValueError: pass
                tmp.append(i)
            self.coords = tuple(tmp)
        else:
            raise TypeError('Atom: Invalid coordinates passed')


    def str(self, symbol=True, suffix='', unit='angstrom', space=14, decimal=6):
        """Return a string representation of this atom.

        Returned string is a single line (no newline characters) that always contains atomic coordinates (and maybe more). Each atomic coordinate is printed using *space* characters, with *decimal* characters reserved for decimal digits. Coordinates values are expressed in *unit*.

        If *symbol* is ``True``, atomic symbol is added at the beginning of the line. If *symbol* is a string, this exact string is printed there.

        *suffix* is an arbitrary string that is appended at the end of returned line. It can contain identifiers in curly brackets (like for example ``f={fragment}``) that will be replaced by values of corresponding attributes (in this case ``self.fragment``). It is done via new string formatting and entire ``self.__dict__`` is passed to formating method. See :ref:`new-string-formatting` for details.

        Example:

            >>> a = Atom(atnum=6, coords=(1,1.5,2))
            >>> print a.str()
                     C      1.000000      1.500000      2.000000
            >>> print a.str(unit='bohr')
                     C      1.889726      2.834589      3.779452
            >>> print a.str(symbol=False)
                  1.000000      1.500000      2.000000
            >>> print a.str(symbol='C2.13')
                 C2.13      1.000000      1.500000      2.000000
            >>> print a.str(suffix='protein1')
                     C      1.000000      1.500000      2.000000 protein1
            >>> a.info = 'membrane'
            >>> print a.str(suffix='subsystem={info}')
                     C      1.000000      1.500000      2.000000 subsystem=membrane

        """
        strformat = '{:>%is}'%space
        numformat = '{:>%i.%if}'%(space,decimal)
        f = lambda x: numformat.format(Units.convert(x, 'angstrom', unit)) if isinstance(x, (int,float)) else strformat.format(str(x))
        if symbol is False:
            return ('{0}{1}{2} '+suffix).format(*map(f,self.coords), **self.__dict__)
        if symbol is True:
            symbol = self.symbol
        return ('{0:>10s}{1}{2}{3} '+suffix).format(symbol, *map(f,self.coords), **self.__dict__)

    def __str__(self):
        """Return a string representation of this atom. Simplified version of :meth:`str` to work as a magic method."""
        return self.str()

    def __iter__(self):
        """Iteration through atom yields coordinates. Thanks to that instances of |Atom| can be passed to any method requiring point or vector as an argument."""
        return iter(self.coords)

    def _setx(self, value): self.coords = (value, self.coords[1], self.coords[2])
    def _sety(self, value): self.coords = (self.coords[0], value, self.coords[2])
    def _setz(self, value): self.coords = (self.coords[0], self.coords[1], value)
    def _getx(self): return self.coords[0]
    def _gety(self): return self.coords[1]
    def _getz(self): return self.coords[2]
    x = property(_getx, _setx)
    y = property(_gety, _sety)
    z = property(_getz, _setz)

    def _getsymbol(self):
        return PT.get_symbol(self.atnum)
    def _setsymbol(self, symbol):
        self.atnum = PT.get_atomic_number(symbol)
    symbol = property(_getsymbol, _setsymbol)

    def _getmass(self):
        return PT.get_mass(self.atnum)
    mass = property(_getmass)

    def _getradius(self):
        return PT.get_radius(self.atnum)
    radius = property(_getradius)

    def _getconnectors(self):
        return PT.get_connectors(self.atnum)
    connectors = property(_getconnectors)

    def translate(self, vector, unit='angstrom'):
        """Move this atom in space by *vector*, expressed in *unit*.

        *vector* should be an iterable container of length 3 (usually tuple, list or numpy array). *unit* describes unit of values stored in *vector*.

        This method requires all coordinates to be numerical values, :exc:`~exceptions.TypeError` is raised otherwise.
        """
        ratio = Units.conversion_ratio(unit, 'angstrom')
        self.coords = tuple(i + j*ratio for i,j in zip(self, vector))


    def move_to(self, point, unit='angstrom'):
        """Move this atom to a given *point* in space, expressed in *unit*.

        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*.

        This method requires all coordinates to be numerical values, :exc:`~exceptions.TypeError` is raised otherwise.
        """
        ratio = Units.conversion_ratio(unit, 'angstrom')
        self.coords = tuple(i*ratio for i in point)


    def distance_to(self, point, unit='angstrom', result_unit='angstrom'):
        """Measure the distance between this atom and *point*.

        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*. Returned value is expressed in *result_unit*.

        This method requires all coordinates to be numerical values, :exc:`~exceptions.TypeError` is raised otherwise.
        """
        ratio = Units.conversion_ratio(unit, 'angstrom')
        res = 0.0
        for i,j in zip(self,point):
            res += (i - j*ratio)**2
        return Units.convert(math.sqrt(res), 'angstrom', result_unit)


    def vector_to(self, point, unit='angstrom', result_unit='angstrom'):
        """Calculate a vector from this atom to *point*.

        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*. Returned value is expressed in *result_unit*.

        This method requires all coordinates to be numerical values, :exc:`~exceptions.TypeError` is raised otherwise.
        """
        ratio = Units.conversion_ratio(unit, 'angstrom')
        resultratio = Units.conversion_ratio('angstrom', result_unit)
        return tuple((i*ratio-j)*resultratio for i,j in zip(point, self))


    def angle(self, point1, point2, point1unit='angstrom', point2unit='angstrom',result_unit='radian'):
        """Calculate an angle between vectors pointing from this atom to *point1* and *point2*.

        *point1* and *point2* should be iterable containers of length 3 (for example: tuple, |Atom|, list, numpy array). Values stored in them are expressed in, respectively, *point1unit* and *point2unit*. Returned value is expressed in *result_unit*.

        This method requires all coordinates to be numerical values, :exc:`~exceptions.TypeError` is raised otherwise.
        """
        num = numpy.dot(self.vector_to(point1, point1unit), self.vector_to(point2, point2unit))
        den = self.distance_to(point1, point1unit) * self.distance_to(point2, point2unit)
        return Units.convert(math.acos(num/den), 'radian', result_unit)

    def rotate(self, matrix):
        """Rotate this atom according to rotation *matrix*.

        *matrix* should be a container with 9 numerical values. It can be a list (tuple, numpy array etc.) listing matrix elements row-wise, either flat (``[1,2,3,4,5,6,7,8,9]``) or in two-level fashion (``[[1,2,3],[4,5,6],[7,8,9]]``).

        .. note::

            This method does not check if supplied matrix is a proper rotation matrix.
        """
        matrix = numpy.array(matrix).reshape(3,3)
        self.coords = tuple(numpy.dot(matrix, numpy.array(self.coords)))



#===================================================================================================
#===================================================================================================
#===================================================================================================



class Bond (object):
    """A class representing a bond between two atoms.

    An instance of this class has the following attributes:
        *   ``atom1`` and ``atom2`` -- two instances of |Atom| that form this bond
        *   ``order`` -- order of the bond. It is either an integer number or the floating point value stored in ``Bond.AR``, indicating aromatic bond
        *   ``mol`` -- a |Molecule| this bond belongs to
        *   ``properties`` -- a |Settings| instance storing all other  information about this bond (initially it is populated with *\*\*other* keyword arguments passed to the constructor)

    .. note::

        Newly created bond is **not** added to ``atom1.bonds`` or ``atom2.bonds``. Storing information about |Bond| in |Atom| is relevant only in the context of the whole |Molecule|, so this information is updated by :meth:`~Molecule.add_bond`.

    """
    AR = 1.5
    def __init__(self, atom1, atom2, order=1, mol=None, **other):
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        self.mol = mol
        self.properties = Settings(other)


    def __str__(self):
        """Return string representation of this bond."""
        return '(%s)--%1.1f--(%s)'%(str(self.atom1), self.order, str(self.atom2))


    def __iter__(self):
        """Iterate over bonded atoms (``atom1`` first, then ``atom2``)."""
        yield self.atom1
        yield self.atom2


    def is_aromatic(self):
        """Check if this bond is aromatic."""
        return self.order == Bond.AR


    def length(self, unit='angstrom'):
        """Return bond's length, expressed in *unit*."""
        return self.atom1.distance_to(self.atom2, result_unit=unit)


    def other_end(self, atom):
        """Return the atom on the other end of this bond with respect to *atom*.

        *atom* has to be either ``atom1`` or ``atom2``, otherwise an exception is raised.
        """
        if atom is self.atom1:
            return self.atom2
        elif atom is self.atom2:
            return self.atom1
        else:
            raise MoleculeError('Bond.other_end: invalid atom passed')


    def resize(self, atom, length, unit='angstrom'):
        """Change the length of the bond to *length*.

        This method works in the following way: one of two atoms forming this bond is moved along the bond in such a way that new length is *length*, in *unit* (direction of the bond in space does not change). Atom indicated by *atom* has to be one of bond's atoms and it is the atom that is **not** moved.
        """
        ratio = 1.0 - Units.convert(length, unit, 'angstrom')/self.length()
        moving = self.other_end(atom)
        moving.translate(tuple(i*ratio for i in moving.vector_to(atom)))



#===================================================================================================
#===================================================================================================
#===================================================================================================



class Molecule (object):
    """A class representing basic molecule object.

    An instance of this class has the following attributes:
        *   ``atoms`` -- a list of |Atom| objects that belong to this molecule
        *   ``bonds`` -- a list of |Bond| objects between atoms listed in ``atoms``
        *   ``lattice`` -- a list of lattice vectors, in case of periodic structures
        *   ``properties`` -- a |Settings| instance storing all other information about this molecule

    .. note::

        Each |Atom| in ``atoms`` list and each |Bond| in ``bonds`` list has a reference to the parent molecule. Moreover, each atom stores the list of bonds it's a part of and each bond stores references to atoms it bonds. That creates a complex net of references between objects that are part of a molecule. Consistency of this data is crucial for proper functioning of many methods. Because of that it is advised not to modify contents of ``atoms`` and ``bonds`` by hand. When you need to alter your molecule, methods :meth:`add_atom`, :meth:`delete_atom`, :meth:`add_bond` and :meth:`delete_bond` can be used to ensure that all these references are updated properly.

    Creating a |Molecule| object for your calculation can be done in two ways. You can start with an empty molecule and manually add all atoms (and bonds, if needed)::

        >>> mol = Molecule()
        >>> mol.add_atom(Atom(atnum=1, coords=(0,0,0)))
        >>> mol.add_atom(Atom(atnum=1, coords=(d,0,0)))

    This approach can be useful for building small molecules, especially if you wish to parametrize some of atomic coordinates (like in :ref:`simple_example`), but in general it's not very practical. Usually one wants to import atomic coordinates from some external file::

        >>> mol = Molecule('xyz/Benzene.xyz')

    Constructor of a |Molecule| object accepts three arguments that can be used to supply this information from a file in your filesystem. *filename* should be a string with a path (absolute or relative) to such a file. *inputformat* describes the format of the file. Currently, the following formats are supported: ``xyz``, ``mol``, ``mol2`` and ``pdb``. If *inputformat* argument is not supplied, PLAMS will try to deduce it by examining the extension of the provided file, so in most of cases it is not needed to use *inputformat*, if only the file has the proper extension. Some formats (``xyz`` and ``pdb``) allow to store more than one geometry of a particular molecule within a single file. In such cases *geometry* argument can be used to indicate which (in order of appearance in the file) geometry to import. *other* keyword arguments passed to the constructor are used to populate ``properties`` |Settings|.

    If a |Molecule| is initialized from an external file, the path to this file (*filename* argument) is stored in ``properties.source``. The base name of the file without extension is kept in ``properties.name``.

    It is also possible to write a molecule to a file in one of the formats mentioned above. See :meth:`write` for details.

    ``lattice`` attribute is used to store information about lattice vectors in case of periodic structures. Some job types (|BANDJob|, |DFTBJob|) will automatically use that data while constructing input files. ``lattice`` should be a list of up to 3 vectors (for different types of periodicity: chain, slab or bulk), each of which needs to be a list or a tuple of 3 numbers.

    Lattice vectors can be directly read and written to ``xyz`` files using the following convention (please mind the fact that this is an unofficial extension to the XYZ format)::

        3

            H      0.000000      0.765440     -0.008360
            O      0.000000      0.000000      0.593720
            H      0.000000     -0.765440     -0.008360
        VEC1       3.000000      0.000000      0.000000
        VEC2       0.000000      3.000000      0.000000
        VEC3       0.000000      0.000000      3.000000

    For 1D (2D) periodicity please supply only ``VEC1`` (``VEC1`` and ``VEC2``). Writing lattice vectors to ``xyz`` files can be disabled by simply reseting the ``lattice`` attribute::

        >>> mol.lattice = []

    |hspace|

    Below the detailed description of available methods is presented. Many of these methods require passing atoms belonging to the molecule as arguments. It can by done by using a reference to an |Atom| object present it ``atoms`` list, but not by passing a number of an atom (its position within ``atoms`` list). Unlike some other tools, PLAMS does not use integer numbers as primary identifiers of atoms. It is done to prevent problems when atoms within a molecule are reordered or some atoms are deleted. References to |Atom| or |Bond| objects can be obtained directly from ``atoms`` or ``bonds`` lists, or with dictionary-like bracket notation::

        >>> mol = Molecule('xyz/Ammonia.xyz')
        >>> mol.guess_bonds()
        >>> print mol
          Atoms:
            1         H      0.942179      0.000000     -0.017370
            2         H     -0.471089      0.815951     -0.017370
            3         N      0.000000      0.000000      0.383210
            4         H     -0.471089     -0.815951     -0.017370
          Bonds:
           (1)--1.0--(3)
           (2)--1.0--(3)
           (3)--1.0--(4)
        >>> at = mol[1]
        >>> print at
                 H      0.942179      0.000000     -0.017370
        >>> b = mol[(1,3)]
        >>> print b
        (         H      0.942179      0.000000     -0.017370 )--1.0--(         N      0.000000      0.000000      0.383210 )
        >>> b = mol[(1,4)]
        >>> print b
        None

    .. note::

        Numbering of atoms within a molecule starts with 1.


    However, if you feel more familiar with identifying atoms by natural numbers, you can use :meth:`set_atoms_id` to equip each atom of the molecule with ``id`` attribute equal to atom's position within ``atoms`` list. This method can also be helpful to track changes in your molecule during tasks that can reorder atoms.
    """


    def __init__(self, filename=None, inputformat=None, geometry=1, **other):
        self.atoms = []
        self.bonds = []
        self.lattice = []
        self.properties = Settings(other)

        if filename is not None :
            self.read(filename, inputformat, geometry)
            self.properties.source = filename
            self.properties.name = os.path.splitext(os.path.basename(filename))[0]


#===================================================================================================
#==== Atoms/bonds manipulation =====================================================================
#===================================================================================================


    def copy(self, atoms=None):
        """Return a copy of this molecule. New molecule has atoms, bonds and all other components distinct from original molecule (it is so called "deep copy").

        By default the entire molecule is copied. It is also possible to copy only some part of the molecule, indicated by *atoms* argument. It should be a list of atoms that **belong to this molecule**. Only these atoms, together with any bonds between them, are copied and included in the returned molecule.
        """
        if atoms is None:
            return copy.deepcopy(self)
        for at in self.atoms:
            at._stay = False
        for at in atoms:
            at._stay = True
        ret = copy.deepcopy(self)
        for at in reversed(ret.atoms):
            if at._stay is False:
                ret.delete_atom(at)
            del at._stay
        for at in self.atoms:
            del at._stay
        return ret


    def add_atom(self, atom, adjacent=None):
        """Add new *atom* to this molecule.

        *atom* should be an |Atom| instance that does not belong to the molecule. Bonds between the new atom and other atoms of the molecule can be automatically added based on *adjacent* argument. It should be a list describing atoms of the molecule that the new atom is connected to. Each element of *adjacent* list can either be a pair ``(Atom, order)`` to indicate new bond's order (use ``Bond.AR`` for aromatic bonds) or an |Atom| instance (a single bond is inserted in this case).

        Example::

            >>> mol = Molecule() #create an empty molecule
            >>> h1 = Atom(symbol='H', coords=(1.0, 0.0, 0.0))
            >>> h2 = Atom(symbol='H', coords=(-1.0, 0.0, 0.0))
            >>> o = Atom(symbol='O', coords=(0.0, 1.0, 0.0))
            >>> mol.add_atom(h1)
            >>> mol.add_atom(h2)
            >>> mol.add_atom(o)
            >>> mol.add_atom(Atom(symbol='C', coords=(0.0, 0.0, 0.0)), adjacent=[h1, h2, (o,2)])

        """
        self.atoms.append(atom)
        atom.mol = self
        if adjacent is not None:
            for adj in adjacent:
                if isinstance(adj, tuple):
                    self.add_bond(atom, adj[0], adj[1])
                else:
                    self.add_bond(atom, adj)


    def delete_atom(self, atom):
        """Delete *atom* from this molecule.

        *atom* should be an |Atom| instance that belongs to the molecule. All bonds containing this atom are removed too.

        Examples::

            >>> #delete all hydrogens
            >>> mol = Molecule('protein.pdb')
            >>> hydrogens = [atom for atom in mol if atom.atnum == 1]
            >>> for i in hydrogens: mol.delete_atom(i)

        ::

            >>> #delete first two atoms
            >>> mol = Molecule('geom.xyz')
            >>> mol.delete_atom(mol[1])
            >>> mol.delete_atom(mol[1]) #since the second atom of original molecule is now the first

        """
        if atom.mol != self:
            raise MoleculeError('delete_atom: passed atom should belong to the molecule')
        try:
            self.atoms.remove(atom)
        except:
            raise MoleculeError('delete_atom: invalid argument passed as atom')
        for b in reversed(atom.bonds):
            self.delete_bond(b)


    def add_bond(self, arg1, arg2=None, order=1):
        """Add new bond to this molecule.

        This method can be used in two different ways. You can call it with just one argument being a |Bond| instance (other arguments are then ignored)::

            >>> b = Bond(mol[2], mol[4], order=Bond.AR) #create aromatic bond between 2nd and 4th atom
            >>> mol.add_bond(b)

        Other way is to pass two atoms (and possibly bond order) and new |Bond| object will be created automatically::

            >>> mol.add_bond(mol[2], mol[4], order=Bond.AR)

        In both cases atoms that are to be bond have to belong to the molecule, otherwise an exception is raised.
        """
        if isinstance(arg1, Atom) and isinstance(arg2, Atom):
            newbond = Bond(arg1, arg2, order=order)
        elif isinstance(arg1, Bond):
            newbond = arg1
        else:
            raise MoleculeError('add_bond: invalid arguments passed')

        if newbond.atom1.mol == self and newbond.atom2.mol == self:
            newbond.mol = self
            self.bonds.append(newbond)
            newbond.atom1.bonds.append(newbond)
            newbond.atom2.bonds.append(newbond)
        else:
            raise MoleculeError('add_bond: bonded atoms have to belong to the molecule')


    def delete_bond(self, arg1, arg2=None):
        """Delete bond from this molecule

        Just like :meth:`add_bond`, this method accepts either a single argument that is a |Bond| instance, or two arguments being instances of |Atom|. In both cases objects used as arguments have to belong to the molecule.
        """
        if isinstance(arg1, Atom) and isinstance(arg2, Atom):
            delbond = self.find_bond(arg1, arg2)
        elif isinstance(arg1, Bond):
            delbond = arg1
        else:
            raise MoleculeError('delete_bond: invalid arguments passed')
        if delbond in self.bonds:
            delbond.mol = None
            self.bonds.remove(delbond)
            delbond.atom1.bonds.remove(delbond)
            delbond.atom2.bonds.remove(delbond)


    def delete_all_bonds(self):
        """Delete all bonds from the molecule."""
        for b in reversed(self.bonds):
            self.delete_bond(b)


    def find_bond(self, atom1, atom2):
        """Find and return a bond between *atom1* and *atom2*. Both atoms have to belong to the molecule. If a bond between chosen atoms does not exist, ``None`` is returned."""
        if atom1.mol != self or atom2.mol != self:
            raise MoleculeError('find_bond: atoms passed as arguments have to belong to the molecule')
        for b in atom1.bonds:
            if atom2 is b.other_end(atom1):
                return b
        return None


    def set_atoms_id(self):
        """Equip each atom of this molecule with ``id`` attribute equal to its position within ``atoms`` list."""
        for i,at in enumerate(self.atoms):
            at.id = i+1

    def unset_atoms_id(self):
        """Delete ``id`` attributes of all atoms."""
        for at in self.atoms:
            try:
                del at.id
            except AttributeError:
                pass


    def neighbors(self, atom):
        """Return a list of neighbors of *atom* within this molecule.

        *atom* has to belong to the molecule. Returned list follows the same order as ``bonds`` list of *atom*.
        """
        if atom.mol != self:
            raise MoleculeError('neighbors: passed atom should belong to the molecule')
        return [b.other_end(atom) for b in atom.bonds]


    def separate(self):
        """Separate this molecule into connected components.

        Returned is a list of new |Molecule| objects (all atoms and bonds are disjoint with original molecule). Each element of this list is identical to one connected component of the base molecule. A connected component is a subset of atoms such that there exists a path (along one or more bonds) between any two atoms.

        Example::

            >>> mol = Molecule('/xyz_dimers/NH3-H2O.xyz')
            >>> mol.guess_bonds()
            >>> print(mol)
              Atoms:
                1         N     -1.395591     -0.021564      0.000037
                2         H     -1.629811      0.961096     -0.106224
                3         H     -1.862767     -0.512544     -0.755974
                4         H     -1.833547     -0.330770      0.862307
                5         O      1.568501      0.105892      0.000005
                6         H      0.606736     -0.033962     -0.000628
                7         H      1.940519     -0.780005      0.000222
              Bonds:
               (5)--1.0--(7)
               (5)--1.0--(6)
               (1)--1.0--(3)
               (1)--1.0--(4)
               (1)--1.0--(2)
            >>> x = mol.separate()
            >>> for i in x: print(i)
              Atoms:
                1         N     -1.395591     -0.021564      0.000037
                2         H     -1.629811      0.961096     -0.106224
                3         H     -1.862767     -0.512544     -0.755974
                4         H     -1.833547     -0.330770      0.862307
              Bonds:
               (1)--1.0--(3)
               (1)--1.0--(4)
               (1)--1.0--(2)

              Atoms:
                1         O      1.568501      0.105892      0.000005
                2         H      0.606736     -0.033962     -0.000628
                3         H      1.940519     -0.780005      0.000222
              Bonds:
               (1)--1.0--(3)
               (1)--1.0--(2)

        """
        frags = []
        clone = self.copy()
        for at in clone:
            at._visited = False

        def dfs(v, mol):
            v._visited = True
            v.mol = mol
            for e in v.bonds:
                e.mol = mol
                u = e.other_end(v)
                if not u._visited:
                    dfs(u, mol)

        for src in clone.atoms:
            if not src._visited:
                m = Molecule()
                dfs(src, m)
                frags.append(m)

        for at in clone.atoms:
            del at._visited
            at.mol.atoms.append(at)
        for b in clone.bonds:
            b.mol.bonds.append(b)

        return frags


    def guess_bonds(self):
        """Try to guess bonds in the molecule based on types and positions of atoms.

        All previously existing bonds are removed. New bonds are generated based on interatomic distances and information about maximal number of bonds for each atom type (``connectors`` property, taken from |PeriodicTable|).

        The problem of finding molecular bonds for a given set of atoms in space does not have a general solution, especially considering the fact the chemical bond is itself not a precisely defined concept. For every method, no matter how sophisticated, there will always be corner cases for which the method produces disputable results. Moreover, depending on the context (area of application) the desired solution for a particular geometry may vary. Please do not treat this method as an oracle always providing proper solution. Algorithm used here gives very good results for geometries that are not very far from optimal geometry, especially consisting of lighter atoms. All kinds of organic molecules, including aromatic ones, usually work very well. Problematic results can emerge for transition metal complexes, transition states, incomplete molecules etc.

        The algorithm used scales as *n log n* where *n* is the number of atoms.

        .. warning::

            This method works reliably only for geometries representing complete molecules. If some atoms are missing (for example, a protein without hydrogens) the resulting set of bonds would usually contain more bonds or bonds with higher order than expected.

        """

        def element(order, ratio, atom1, atom2):
            eford = order
            if order == 1.5:
                eford = 1.15
            elif order == 1 and {atom1.symbol, atom2.symbol} == {'C', 'N'}:
                eford = 1.11
            return ((eford+0.9)*ratio, order, ratio, atom1, atom2)

        self.delete_all_bonds()

        dmax = 1.28

        cubesize = dmax*2.1*max([at.radius for at in self.atoms])

        cubes = {}
        for i,at in enumerate(self.atoms):
            at._id = i+1
            at.free = at.connectors
            at.cube = tuple(map(lambda x: int(math.floor(x/cubesize)), at.coords))
            if at.cube in cubes:
                cubes[at.cube].append(at)
            else:
                cubes[at.cube] = [at]

        neighbors = {}
        for cube in cubes:
            neighbors[cube] = []
            for i in range(cube[0]-1, cube[0]+2):
                for j in range(cube[1]-1, cube[1]+2):
                    for k in range(cube[2]-1, cube[2]+2):
                        if (i,j,k) in cubes:
                            neighbors[cube] += cubes[(i,j,k)]

        heap = []
        for at1 in self.atoms:
            if at1.free > 0:
                for at2 in neighbors[at1.cube]:
                    if (at2.free > 0) and (at1._id < at2._id):
                        ratio = at1.distance_to(at2)/(at1.radius+at2.radius)
                        if (ratio < dmax):
                            heap.append(element(0, ratio, at1, at2))
                            #I hate to do this, but I guess there's no other way :/ [MH]
                            if (at1.atnum == 16 and at2.atnum == 8):
                                at1.free = 6
                            elif (at2.atnum == 16 and at1.atnum == 8):
                                at2.free = 6
                            elif (at1.atnum == 7):
                                at1.free += 1
                            elif (at2.atnum == 7):
                                at2.free += 1
        heapq.heapify(heap)

        for at in self.atoms:
            if at.atnum == 7:
                if at.free > 6:
                    at.free = 4
                else:
                    at.free = 3

        while heap:
            val, o, r, at1, at2 = heapq.heappop(heap)
            step = 1 if o in [0,2] else 0.5
            if at1.free >= step and at2.free >= step:
                o += step
                at1.free -= step
                at2.free -= step
                if o < 3:
                    heapq.heappush(heap, element(o,r,at1,at2))
                else:
                    self.add_bond(at1,at2,o)
            elif o > 0:
                if o == 1.5:
                    o = Bond.AR
                self.add_bond(at1,at2,o)

        def dfs(atom, par):
            atom.arom += 1000
            for b in atom.bonds:
                oe = b.other_end(atom)
                if b.is_aromatic() and oe.arom < 1000:
                    if oe.arom > 2:
                        return False
                    if par and oe.arom == 1:
                        b.order = 2
                        return True
                    if dfs(oe, 1-par):
                        b.order = 1 + par
                        return True

        for at in self.atoms:
            at.arom = len(list(filter(Bond.is_aromatic, at.bonds)))

        for at in self.atoms:
            if at.arom == 1:
                dfs(at, 1)

        for at in self.atoms:
            del at.cube,at.free,at._id,at.arom




#===================================================================================================
#==== Geometry operations ==========================================================================
#===================================================================================================



    def translate(self, vector, unit='angstrom'):
        """Move this molecule in space by *vector*, expressed in *unit*.

        *vector* should be an iterable container of length 3 (usually tuple, list or numpy array). *unit* describes unit of values stored in *vector*.
        """
        for at in self.atoms:
            at.translate(vector, unit)


    def rotate(self, matrix):
        """Rotate this molecule according to rotation *matrix*.

        *matrix* should be a container with 9 numerical values. It can be a list (tuple, numpy array etc.) listing matrix elements row-wise, either flat (``[1,2,3,4,5,6,7,8,9]``) or in two-level fashion (``[[1,2,3],[4,5,6],[7,8,9]]``).

        .. note::

            This method does not check if supplied matrix is a proper rotation matrix.
        """
        for at in self.atoms:
            at.rotate(matrix)


    def rotate_bond(self, bond, atom, angle, unit='radian'):
        """Rotate given *bond* by an *angle* expressed in *unit*.

        *bond* should be chosen in such a way, that it divides the molecule into two parts (using a bond being part of a ring results in an error). *atom* has to belong to *bond* and is used to pick which "half" of the molecule is rotated. Positive angle denotes counterclockwise rotation (looking along the bond, from the stationary part of the molecule).
        """
        if atom not in bond:
            raise MoleculeError('rotate_bond: atom has to belong to the bond')

        atoms_to_rotate = {atom}

        def dfs(v):
            for e in v.bonds:
                if e is not bond:
                    u = e.other_end(v)
                    if u not in atoms_to_rotate:
                        atoms_to_rotate.add(u)
                        dfs(u)

        dfs(atom)

        if len(atoms_to_rotate) == len(self):
            raise MoleculeError('rotate_bond: chosen bond does not divide molecule')

        other_end = bond.other_end(atom)
        v = numpy.array(other_end.vector_to(atom))
        v /= numpy.linalg.norm(v)

        W = numpy.array([[0, -v[2], v[1]],
                         [v[2], 0, -v[0]],
                         [-v[1], v[0], 0]])

        angle = Units.convert(angle, unit, 'radian')
        a1 = math.sin(angle)
        a2 = 2 * math.pow(math.sin(0.5 * angle), 2)

        rotmat = numpy.identity(3) + a1 * W + a2 * numpy.dot(W,W)

        trans = numpy.array(other_end.vector_to((0,0,0)))
        for at in atoms_to_rotate:
            at.translate(trans)
            at.rotate(rotmat)
            at.translate(-trans)


    def closest_atom(self, point, unit='angstrom'):
        """Return the atom of this molecule that is the closest one to some *point* in space.

        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*.
        """
        dist = float('inf')
        for at in self.atoms:
            newdist = at.distance_to(point, unit=unit)
            if newdist < dist:
                dist = newdist
                ret = at
        return ret


    def distance_to_point(self, point, unit='angstrom', result_unit='angstrom'):
        """Calculate the distance between this molecule and some *point* in space (distance between *point* and :meth:`closest_atom`).

        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*. Returned value is expressed in *result_unit*.
        """
        at = self.closest_atom(point, unit)
        return at.distance_to(point, unit, result_uni)


    def distance_to_mol(self, other, result_unit='angstrom', return_atoms=False):
        """Calculate the distance between this molecule and some *other* molecule.

        The distance is measured as the smallest distance between a pair of atoms, one belonging to each of the molecules. Returned distance is expressed in *result_unit*.

        If *return_atoms* is ``False``, only a single number is returned.  If *return_atoms* is ``True``, this method returns a tuple ``(distance, atom1, atom2)`` where ``atom1`` and ``atom2`` are atoms fulfilling the minimal distance, with atom1 belonging to this molecule and atom2 to *other*.
        """
        dist = float('inf')
        for at1 in self.atoms:
            for at2 in other.atoms:
                newdist = (at1.x-at2.x)**2 + (at1.y-at2.y)**2 + (at1.z-at2.z)**2
                if newdist < dist:
                    dist = newdist
                    atom1 = at1
                    atom2 = at2
        res = Units.convert(math.sqrt(dist), 'angstrom', result_unit)
        if return_atoms:
            return res, atom1, atom2
        return res


    def wrap(self, length, angle=2*math.pi, length_unit='angstrom', angle_unit='radian'):
        """wrap(self, length, angle=2*pi, length_unit='angstrom', angle_unit='radian')

        Transform the molecule wrapping its x-axis around z-axis. This method is useful for building nanotubes or molecular wedding rings.

        Atomic coordinates are transformed in the following way:
            *   zzzzz coordinates remain untouched
            *   x axis gets wrapped around the circle centered in the origin of new coordinate system. Each segment of x axis of length *length* ends up as an arc of a circle subtended by an angle *angle*. The radius of this circle is R = *length*/*angle*.
            *   part of the plane between the x axis and the line y=R is transformed into the interior of the circle, with line y=R being squashed into a single point - the center of the circle.
            *   part of the plane above line y=R is dropped
            *   part of the plane below x axis is transformed into outside of the circle
            *   transformation is done in such a way that distances along y axis are preserved

        Before:

        .. image:: _static/wrap.*

        After:

        .. image:: _static/wrap2.*

        """
        length = Units.convert(length, length_unit, 'angstrom')
        angle = Units.convert(angle, angle_unit, 'radian')

        xs = [atom.x for atom in self.atoms]
        if max(xs)-min(xs) > length:
            raise MoleculeError('wrap: x-extension of the molecule is larger than length')

        if angle < 0 or angle > 2*math.pi:
            raise MoleculeError('wrap: angle must be between 0 and 2*pi')

        R = length / angle

        def map_ring(x,y):
            return ((R-y) * math.cos(x/R), (R-y) * math.sin(x/R))

        for at in self.atoms:
            at.x, at.y = map_ring(at.x, at.y)


    def get_center_of_mass(self, unit='angstrom'):
        """Return the center of mass of this molecule (as a tuple). Returned coordinates are expressed in *unit*."""
        center = [0.0,0.0,0.0]
        total_mass = 0.0
        for at in self.atoms:
            total_mass += at.mass
            for i in range(3):
                center[i] += at.mass*at.coords[i]
        for i in range(3):
            center[i] = Units.convert(center[i]/total_mass, 'angstrom', unit)
        return tuple(center)


    def get_mass(self):
        """Return mass of the molecule, expressed in atomic units."""
        return sum([at.mass for at in self.atoms])


    def get_formula(self):
        """Calculate the molecular formula for this molecule.

        Returned value is a single string. It contains simple molecular formula (it only includes atom types and total number of atoms of each type)."""
        atnums = [at.atnum for at in self.atoms]
        s = set(atnums)
        formula = ''
        for i in s:
            formula += PT.get_symbol(i) + str(atnums.count(i))
        return formula



#===================================================================================================
#==== Magic methods ================================================================================
#===================================================================================================



    def __len__(self):
        """Length of a molecule is the number of atoms."""
        return len(self.atoms)

    def __str__(self):
        """Return string representation of this molecule.

        Information about atoms are printed in ``xyz`` format fashion -- each atom in a separate, enumerated line. Then, if the molecule contains any bonds, they are printed. Each bond is printed in a separate line, with information about both atoms and bond order. Example::

                  Atoms:
                    1         N       0.00000       0.00000       0.38321
                    2         H       0.94218       0.00000      -0.01737
                    3         H      -0.47109       0.81595      -0.01737
                    4         H      -0.47109      -0.81595      -0.01737
                  Bonds:
                    (1)----1----(2)
                    (1)----1----(3)
                    (1)----1----(4)
        """
        s = '  Atoms: \n'
        for i,atom in enumerate(self.atoms):
            s += ('%5i'%(i+1)) + str(atom) + '\n'
        if len(self.bonds) > 0:
            for j,atom in enumerate(self.atoms):
                atom._tmpid = j+1
            s += '  Bonds: \n'
            for bond in self.bonds:
                s += '   (%d)--%1.1f--(%d)\n'%(bond.atom1._tmpid, bond.order, bond.atom2._tmpid)
            for atom in self.atoms:
                del atom._tmpid
        if self.lattice:
            s += '  Lattice:\n'
            for vec in self.lattice:
               s += '    %10.6f %10.6f %10.6f\n'%vec
        return s


    def __iter__(self):
        """Iterate over atoms."""
        return iter(self.atoms)


    def __getitem__(self, key):
        """Bracket notation can be used to access atoms or bonds directly.

        If *key* is a single int (``mymol[i]``), return i-th atom of this molecule. If *key* is a pair of ints (``mymol[(i,j)]``), return bond between i-th and j-th atom (``None`` if such a bond does not exist).

        This method is read only (things like ``mymol[3] = Atom(...)`` are forbidden). Numbering of atoms withing a molecule starts with 1.

        """
        if isinstance(key, int):
            if key == 0:
                raise MoleculeError('Numbering of atoms starts with 1')
            if key < 0:
                return self.atoms[key]
            return self.atoms[key-1]
        if isinstance(key, tuple) and len(key) == 2:
            if key[0] == 0 or key[1] == 0:
                raise MoleculeError('Numbering of atoms starts with 1')
            return self.find_bond(self.atoms[key[0]-1], self.atoms[key[1]-1])

    def __add__(self, other):
        """Create a new molecule that is a sum of this molecule and *other*::

        >>> newmol = mol1 + mol2

        The new molecule has atoms, bonds and all other elements distinct from both components. ``properties`` of ``newmol`` are ``properties`` of ``mol1`` :meth:`soft_updated<scm.plams.settings.Settings.soft_update>` with ``properties`` of ``mol2``.
        """
        m = self.copy()
        m += other
        return m


    def __iadd__(self, other):
        """Add *other* molecule to this one::

        >>> protein += water

        All atoms and bonds present in *other* are copied and copies are added to this molecule. ``properties`` of this molecule are :meth:`soft_updated<scm.plams.settings.Settings.soft_update>` with ``properties`` of *other*.
        """
        othercopy = other.copy()
        self.atoms += othercopy.atoms
        self.bonds += othercopy.bonds
        for atom in self.atoms:
            atom.mol = self
        for bond in self.bonds:
            bond.mol = self
        self.properties.soft_update(othercopy.properties)
        return self


    def __copy__(self):
        return self.copy()



#===================================================================================================
#==== File/format IO ===============================================================================
#===================================================================================================



    def readxyz(self, f, frame):

        def newatom(line):
            lst = line.split()
            shift = 1 if (len(lst) > 4 and lst[0] == str(i)) else 0
            num = lst[0+shift]
            if isinstance(num, str):
                num = PT.get_atomic_number(num)
            self.add_atom(Atom(atnum=num, coords=(lst[1+shift],lst[2+shift],lst[3+shift])))

        def newlatticevec(line):
            lst = line.split()
            self.lattice.append((float(lst[1]),float(lst[2]),float(lst[3])))

        fr = frame
        begin, first, nohead = True, True, False
        for line in f:
            if first:
                if line.strip() == '' : continue
                first = False
                try:
                    n = int(line.strip())
                    fr -= 1
                except ValueError:
                    nohead = True
                    newatom(line)
            elif nohead:
                if line.strip() == '' : break
                if 'VEC' in line.upper():
                    newlatticevec(line)
                else:
                    newatom(line)
            elif fr != 0:
                try:
                    n = int(line.strip())
                    fr -= 1
                except ValueError:
                    continue
            else:
                if begin:
                    begin = False
                    i = 1
                    if line:
                        self.properties['comment'] = line.rstrip()
                else:
                    if i <= n:
                        newatom(line)
                        i += 1
                    elif 'VEC' in line.upper():
                       newlatticevec(line)
                    else:
                        break
        if not nohead and fr > 0:
            raise FileError('readxyz: There are only %i frames in %s' % (frame - fr, f.name))


    def writexyz(self, f):
        f.write(str(len(self)) + '\n')
        if 'comment' in self.properties:
            comment = self.properties['comment']
            if isinstance(comment, list):
                comment = comment[0]
            f.write(comment)
        f.write('\n')
        for at in self.atoms:
            f.write(str(at) + '\n')
        for i,vec in enumerate(self.lattice):
            f.write('VEC'+str(i+1) + '%14.6f %14.6f %14.6f\n'%tuple(vec))


    def readmol(self, f, frame):
        if frame != 1:
            raise FileError('readmol: .mol files do not support multiple geometries')

        comment = []
        for i in range(4):
            line = f.readline().rstrip()
            if line:
                spl = line.split()
                if spl[len(spl)-1] == 'V2000':
                    natom = int(spl[0])
                    nbond = int(spl[1])
                    for j in range(natom):
                        atomline = f.readline().split()
                        crd = tuple(map(float, atomline[0:3]))
                        symb = atomline[3]
                        try:
                            num = PT.get_atomic_number(symb)
                        except PTError:
                            num = 0
                        self.add_atom(Atom(atnum=num, coords=crd))
                    for j in range(nbond):
                        bondline = f.readline().split()
                        at1 = self.atoms[int(bondline[0]) - 1]
                        at2 = self.atoms[int(bondline[1]) - 1]
                        ordr = int(bondline[2])
                        if ordr == 4:
                            ordr = Bond.AR
                        self.add_bond(Bond(atom1=at1, atom2=at2, order=ordr))
                    break
                elif spl[len(spl)-1] == 'V3000':
                    raise FileError('readmol: Molfile V3000 not supported. Please convert')
                else:
                    comment.append(line)
        if comment:
            self.properties['comment'] = comment



    def writemol(self, f):
        commentblock = ['\n']*3
        if 'comment' in self.properties:
            comment = self.properties['comment']
            if isinstance(comment, str):
                commentblock[0] = comment + '\n'
            elif isinstance(comment, list):
                comment = comment[0:3]
                while len(comment) < 3:
                    comment.append('')
                commentblock = [a+b for a,b in zip(comment,commentblock)]
        f.writelines(commentblock)

        self.set_atoms_id()

        f.write('%3i%3i  0  0  0  0  0  0  0  0999 V2000\n' % (len(self.atoms),len(self.bonds)))
        for at in self.atoms:
            f.write('%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0\n' % (at.x,at.y,at.z,at.symbol))
        for bo in self.bonds:
            order = bo.order
            if order == Bond.AR:
                order = 4
            f.write('%3i%3i%3i  0  0  0\n' % (bo.atom1.id,bo.atom2.id,order))
        self.unset_atoms_id()
        f.write('M  END\n')



    def readmol2(self, f, frame):
        if frame != 1:
            raise MoleculeError('readmol: .mol2 files do not support multiple geometries')

        bondorders = {'1':1, '2':2, '3':3, 'am':1, 'ar':Bond.AR, 'du':0, 'un':1, 'nc':0}
        mode = ('', 0)
        for i, line in enumerate(f):
            line = line.rstrip()
            if not line:
                continue
            elif line[0] == '#':
                continue
            elif line[0] == '@':
                line = line.partition('>')[2]
                if not line:
                    raise FileError('readmol2: Error in %s line %i: invalid @ record' % (f.name, str(i+1)))
                mode = (line, i)

            elif mode[0] == 'MOLECULE':
                pos = i - mode[1]
                if pos == 1:
                    self.properties['name'] = line
                elif pos == 3:
                    self.properties['type'] = line
                elif pos == 4:
                    self.properties['charge_type'] = line
                elif pos == 5:
                    self.properties['flags'] = line
                elif pos == 6:
                    self.properties['comment'] = line

            elif mode[0] == 'ATOM':
                spl = line.split()
                if len(spl) < 6:
                    raise FileError('readmol2: Error in %s line %i: not enough values in line' % (f.name, str(i+1)))
                symb = spl[5].partition('.')[0]
                try:
                    num = PT.get_atomic_number(symb)
                except PTError:
                    num = 0
                crd = tuple(map(float, spl[2:5]))
                newatom = Atom(atnum=num, coords=crd, name=spl[1], type=spl[5])
                if len(spl) > 6:
                    newatom.properties['subst_id'] = spl[6]
                if len(spl) > 7:
                    newatom.properties['subst_name'] = spl[7]
                if len(spl) > 8:
                    newatom.properties['charge'] = float(spl[8])
                if len(spl) > 9:
                    newatom.properties['flags'] = spl[9]
                self.add_atom(newatom)

            elif mode[0] == 'BOND':
                spl = line.split()
                if len(spl) < 4:
                    raise FileError('readmol2: Error in %s line %i: not enough values in line' % (f.name, str(i+1)))
                try:
                    atom1 = self.atoms[int(spl[1])-1]
                    atom2 = self.atoms[int(spl[2])-1]
                except IndexError:
                    raise FileError('readmol2: Error in %s line %i: wrong atom ID' % (f.name, str(i+1)))
                newbond = Bond(atom1, atom2, order=bondorders[spl[3]])
                if len(spl) > 4:
                    for flag in spl[4].split('|'):
                        newbond.properties[flag] = True
                self.add_bond(newbond)



    def writemol2(self, f):
        bondorders = ['1','2','3','ar']

        def write_prop(name, obj, separator, space=0, replacement=None):
            form_str = '%-' + str(space) + 's'
            if name in obj.properties:
                f.write(form_str % str(obj.properties[name]))
            elif replacement is not None:
                f.write(form_str % str(replacement))
            f.write(separator)

        f.write('@<TRIPOS>MOLECULE\n')
        write_prop('name', self, '\n')
        f.write('%i %i\n' % (len(self.atoms),len(self.bonds)))
        write_prop('type', self, '\n')
        write_prop('charge_type', self, '\n')
        write_prop('flags', self, '\n')
        write_prop('comment', self, '\n')

        f.write('\n@<TRIPOS>ATOM\n')
        for i,at in enumerate(self.atoms):
            f.write('%5i ' % (i+1))
            write_prop('name', at, ' ', 5, at.symbol+str(i+1))
            f.write('%10.4f %10.4f %10.4f ' % at.coords)
            write_prop('type', at, ' ', 5, at.symbol)
            write_prop('subst_id', at, ' ', 5)
            write_prop('subst_name', at, ' ', 7)
            write_prop('charge', at, ' ', 6)
            write_prop('flags', at, '\n')
            at.id = i+1

        f.write('\n@<TRIPOS>BOND\n')
        for i,bo in enumerate(self.bonds):
            f.write('%5i %5i %5i %4s' % (i+1, bo.atom1.id, bo.atom2.id, bondorders[bo.order]))
            write_prop('flags', bo, '\n')

        self.unset_atoms_id()



    def readpdb(self, f, frame):
        pdb = PDBHandler(f)
        models = pdb.get_models()
        if frame > len(models):
            raise FileError('readpdb: There are only %i frames in %s' % (len(models), f.name))

        symbol_columns = [70,6,7,8]
        for i in models[frame-1]:
            if i.name in ['ATOM  ','HETATM']:
                x = float(i.value[0][24:32])
                y = float(i.value[0][32:40])
                z = float(i.value[0][40:48])
                for n in symbol_columns:
                    symbol = i.value[0][n:n+2].strip()
                    try:
                        atnum = PT.get_atomic_number(symbol)
                        break
                    except PTError:
                        if n == symbol_columns[-1]:
                            raise FileError('readpdb: Unable to deduce the atomic symbol in the following line:\n%s'%(i.name+i.value[0]))
                self.add_atom(Atom(atnum=atnum,coords=(x,y,z)))

        return pdb



    def writepdb(self, f):
        pdb = PDBHandler()
        pdb.add_record(PDBRecord('HEADER'))
        model = []
        for i,at in enumerate(self.atoms):
            s = 'ATOM  %5i                   %8.3f%8.3f%8.3f                      %2s  ' % (i+1,at.x,at.y,at.z,at.symbol.upper())
            model.append(PDBRecord(s))
        pdb.add_model(model)
        pdb.add_record(pdb.calc_master())
        pdb.add_record(PDBRecord('END'))
        pdb.write(f)


    def read(self, filename, inputformat=None, frame=1):
        """Read molecular coordinates from file.

        *filename* should be a string with a path to the file. If *inputformat* is not ``None``, it should be one of supported formats (keys occurring in class attribute ``_readformat``). Otherwise, format of the file is deduced from file's extension (for files without extension `xyz` format is assumed).

        If chosen format allows multiple geometries in a single file, *frame* can be used to pick one of them.
        """
        if inputformat is None:
            fsplit = filename.rsplit('.',1)
            if len(fsplit) == 2:
                inputformat = fsplit[1]
            else:
                inputformat = 'xyz'
        if inputformat in self.__class__._readformat:
            with open(filename, 'rU') as f:
                ret = self._readformat[inputformat](self, f, frame)
            return ret
        else:
            raise MoleculeError('read: Unsupported file format')



    def write(self, filename, outputformat=None):
        """Write molecular coordinates to a file.

        *filename* should be a string with a path to the file. If *outputformat* is not ``None``, it should be one of supported formats (keys occurring in class attribute ``_writeformat``). Otherwise, format of the file is deduced from file's extension (for files without extension `xyz` format is assumed).
        """
        if outputformat is None:
            fsplit = filename.rsplit('.',1)
            if len(fsplit) == 2:
                outputformat = fsplit[1]
            else:
                outputformat = 'xyz'
        if outputformat in self.__class__._writeformat:
            with open(filename, 'w') as f:
                self._writeformat[outputformat](self, f)
        else:
            raise MoleculeError('write: Unsupported file format')

    _readformat = {'xyz':readxyz, 'mol':readmol, 'mol2':readmol2, 'pdb':readpdb}
    _writeformat = {'xyz':writexyz, 'mol':writemol, 'mol2':writemol2, 'pdb': writepdb}

#===================================================================================================
#==== JSON IO ======================================================================================
#===================================================================================================


    def as_dict(self):
        """
        The Molecule information is stored in a dict based on the mol
        `file format <http://onlinelibrarystatic.wiley.com/marvin/help/FF/19693841.html>`
        :returns: JSON object

        """
        # def create_atom_ids(mol):
        #     """
        #     Generate unique identifier for each atom
        #     :parameter mol: Molecule object containing the atoms
        #     :type      mol: |Molecule|
        #     """
        #     ds = []
        #     for i,self.assertTrue() in enumerate(mol.atoms):
        #         idx  = "atom_" + str(i)
        #         ds.append(idx)
        #     return ds

        def create_atom_block(mol):
            """
            In this block are stored the atomic number, coordinates,
            atomic symbols and others properties that are passed as a *Setting* object.
            :parameter mol: Molecule object containing the atoms
            :type      mol: |Molecule|
            :parameter ids: List of unique identifier
            :type      ids: [Int]
            """
            return [{'coords': at.coords, 'symbol': at.symbol,
                     'atnum': at.atnum, 'properties': at.properties.as_dict()} for at in mol.atoms]

        def create_bond_block(mol):
            """
            :parameter mol: Molecule object containing the atoms
            :type      mol: |Molecule|
            :parameter ids: List of unique identifier
            :type      ids: [Int]

            The data list on table bellow is used to store
            information related to the bonds.

            =================   =====================
            Meaning             Value
            ==================  =====================
            First atom number   Int
            Second atom number  Int
            Bond type           Int (float for aromatic)
                                #. Single
                                #. Double
                                #. Triple
                                #. Aromatic
            Properties          Settings
            ==================  =====================
            """
            def get_bond_order(bond):
                if bond.order:
                    return bond.order
                elif bond.AR:
                    return bond.AR
                else:
                    msg = "bond does not contain order attribute"
                    raise AttributeError(msg)

            # Represent the atoms involved in the bond like references
            # to the unique identifiers stored in the ``atomBlock``
            # dict_atom2id = {at: idx for (at, idx) in zip(mol.atoms, ids)}
            bonds_list = []
            for i, b in enumerate(mol.bonds):
                atom1 = mol.atoms.index(b.atom1)
                atom2 = mol.atoms.index(b.atom2)
                order = get_bond_order(b)
                bond  = {'atom1': atom1, 'atom2': atom2,
                         'order': order, 'properties': b.properties.as_dict()}
                bonds_list.append(bond)

            return bonds_list

        d              = dict()
        d["atomBlock"] = create_atom_block(self)
        d["bondBlock"] = create_bond_block(self)
        d["properties"] = self.properties

        return d

    @classmethod
    def from_dict(cls, atomBlock, bondBlock, properties):
        """
        Generate a new Molecule instance using the data stored
        in the dictionary representing the JSON serialized data
        :parameter ds: Dict containing the JSON serialized molecule
        :type      ds: Dict
        :returns: |Molecule|
        """
        # New Molecule instance
        mol     = cls()
        # dict from unique atom identifiers to numeration
        # inside the molecule

        # Reconstruct the Atom instances using the Json data
        for at in atomBlock:
            atnum     = at["atnum"]
            coords    = at["coords"]
            symbol    = at["symbol"]
            props     = at["properties"]
            mol.add_atom(Atom(atnum=atnum, coords=coords, symbol=symbol,
                              **props))
        # Reconstruct the bonds using the internal numeration of the molecule
        # build in  the previous step.
        for b in bondBlock:
            id_atom1 = b["atom1"]
            id_atom2 = b["atom2"]
            atom1    = mol.atoms[id_atom1]
            atom2    = mol.atoms[id_atom2]
            bond = Bond(atom1=atom1, atom2=atom2, order=b["order"],
                        **b["properties"])
            mol.add_bond(bond)

        mol.properties = properties

        return mol
