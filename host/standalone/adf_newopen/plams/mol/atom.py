import math
import numpy as np

from ..core.settings import Settings
from ..tools.periodic_table import PT
from ..tools.units import Units


__all__ = ['Atom']


class Atom:
    """A class representing a single atom in three dimensional space.

    An instance of this class has the following attributes:

    *   ``atnum`` -- atomic number (zero for "dummy atoms")
    *   ``coords`` -- tuple of length 3 storing spatial coordinates
    *   ``bonds`` -- list of bonds (see |Bond|) this atom is a part of
    *   ``mol`` -- |Molecule| this atom belongs to
    *   ``properties`` -- |Settings| instance storing all other information about this atom (initially it is populated with *\*\*other*)

    The above attributes can be accessed either directly or using one of the following properties:

    *   ``x``, ``y``, ``z`` -- allow to read or modify each coordinate separately
    *   ``symbol`` -- allows to read or write the atomic symbol directly. Atomic symbol is not stored as an attribute, instead of that the atomic number (``atnum``) indicates the type of atom. In fact, ``symbol`` this is just a wrapper around ``atnum`` that uses |PeriodicTable| as a translator::

            >>> a = Atom(atnum=8)
            >>> print(a.symbol)
            O
            >>> a.symbol = 'Ca'
            >>> print(a.atnum)
            20

    *   ``mass`` -- atomic mass, obtained from |PeriodicTable|, read only
    *   ``radius`` -- atomic radius, obtained from |PeriodicTable|, read only
    *   ``connectors`` -- number of connectors, obtained from |PeriodicTable|, read only

    .. note::

        When creating a new atom, its type can be chosen either by setting an atomic number or a symbol (``atnum`` and ``symbol`` constructor arguments). The ``symbol`` argument takes precedence -- if it is supplied, the ``atnum`` argument is ignored.

    Values stored in ``coords`` tuple do not necessarily have to be numeric, you can also store any string there. This might come handy for programs that allow parametrization of coordinates in the input file (to enforce some geometry constraints for example)::

            >>> a = Atom(symbol='C', coords=(1,2,3))
            >>> print(a)
                     C       1.00000       2.00000       3.00000
            >>> a.y = 'param1'
            >>> print(a)
                     C       1.00000        param1       3.00000

    However, non-numerical coordinates cannot be used together with some methods (for example :meth:`distance_to` or :meth:`translate`). An attempt to do this raises an exception.

    Internally, atomic coordinates are always expressed in angstroms. Most of methods that read or modify atomic coordinates accept a keyword argument ``unit`` allowing to choose unit in which results and/or arguments are expressed (see |Units| for details). Throughout the entire code angstrom is the default length unit. If you don't specify ``unit`` parameter in any place of your script, all the automatic unit handling described above boils down to occasional multiplication/division by 1.0.
    """
    def __init__(self, atnum=0, symbol=None, coords=None, unit='angstrom', bonds=None, mol=None, **other):
        if symbol is not None:
            self.symbol = symbol
        else:
            self.atnum = atnum
            if atnum == 0: self._dummysymbol = 'Xx'

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


    def str(self, symbol=True, suffix='', suffix_dict={}, unit='angstrom', space=14, decimal=6):
        """Return a string representation of this atom.

        Returned string is a single line (no newline characters) that always contains atomic coordinates (and maybe more). Each atomic coordinate is printed using *space* characters, with *decimal* characters reserved for decimal digits. Coordinates values are expressed in *unit*.

        If *symbol* is ``True``, atomic symbol is added at the beginning of the line. If *symbol* is a string, this exact string is printed there.

        *suffix* is an arbitrary string that is appended at the end of returned line. It can contain identifiers in curly brackets (like for example ``f={fragment}``) that will be replaced by values of corresponding keys from *suffix_dict* dictionary. See :ref:`formatstrings` for details.

        Example:

            >>> a = Atom(atnum=6, coords=(1,1.5,2))
            >>> print(a.str())
                     C      1.000000      1.500000      2.000000
            >>> print(a.str(unit='bohr'))
                     C      1.889726      2.834589      3.779452
            >>> print(a.str(symbol=False))
                  1.000000      1.500000      2.000000
            >>> print(a.str(symbol='C2.13'))
                 C2.13      1.000000      1.500000      2.000000
            >>> print(a.str(suffix='protein1'))
                     C      1.000000      1.500000      2.000000 protein1
            >>> a.properties.info = 'membrane'
            >>> print(a.str(suffix='subsystem={info}', suffix_dict=a.properties))
                     C      1.000000      1.500000      2.000000 subsystem=membrane

        """
        strformat = '{:>%is}'%space
        numformat = '{:>%i.%if}'%(space,decimal)
        f = lambda x: numformat.format(Units.convert(x, 'angstrom', unit)) if isinstance(x, (int,float)) else strformat.format(str(x))
        if symbol is False:
            return ('{0} {1} {2} '+suffix).format(*map(f,self.coords), **suffix_dict)
        if symbol is True:
            symbol = self.symbol
        return ('{0:>10s}{1} {2} {3} '+suffix).format(symbol, *map(f,self.coords), **suffix_dict)


    def __str__(self):
        """Return a string representation of this atom. Simplified version of :meth:`str` to work as a magic method."""
        return self.str()


    def __iter__(self):
        """Iteration through atom yields coordinates. Thanks to that instances of |Atom| can be passed to any method requiring as an argument a point or a vector in 3D space."""
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
        if self.atnum == 0:
            return self._dummysymbol
        else:
            return PT.get_symbol(self.atnum)
    def _setsymbol(self, symbol):
        if symbol.lower().capitalize() in PT.dummysymbols:
            self.atnum = 0
            self._dummysymbol = symbol.lower().capitalize()
        else:
            self.atnum = PT.get_atomic_number(symbol)
            self._dummysymbol = None
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

    def _ismetallic(self):
        return PT.get_metallic(self.atnum)
    is_metallic = property(_ismetallic)

    def _iselectronegative(self):
        return PT.get_electronegative(self.atnum)
    is_electronegative = property(_iselectronegative)



    def translate(self, vector, unit='angstrom'):
        """Move this atom in space by *vector*, expressed in *unit*.

        *vector* should be an iterable container of length 3 (usually tuple, list or numpy array). *unit* describes unit of values stored in *vector*.

        This method requires all atomic coordinates to be numerical values, :exc:`~exceptions.TypeError` is raised otherwise.
        """
        ratio = Units.conversion_ratio(unit, 'angstrom')
        self.coords = tuple(i + j*ratio for i,j in zip(self, vector))


    def move_to(self, point, unit='angstrom'):
        """Move this atom to a given *point* in space, expressed in *unit*.

        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*.

        This method requires all atomic coordinates to be numerical values, :exc:`~exceptions.TypeError` is raised otherwise.
        """
        ratio = Units.conversion_ratio(unit, 'angstrom')
        self.coords = tuple(i*ratio for i in point)


    def distance_to(self, point, unit='angstrom', result_unit='angstrom'):
        """Measure the distance between this atom and *point*.

        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*. Returned value is expressed in *result_unit*.

        This method requires all atomic coordinates to be numerical values, :exc:`~exceptions.TypeError` is raised otherwise.
        """
        ratio = Units.conversion_ratio(unit, 'angstrom')
        res = 0.0
        for i,j in zip(self,point):
            res += (i - j*ratio)**2
        return Units.convert(math.sqrt(res), 'angstrom', result_unit)


    def vector_to(self, point, unit='angstrom', result_unit='angstrom'):
        """Calculate a vector from this atom to *point*.

        *point* should be an iterable container of length 3 (for example: tuple, |Atom|, list, numpy array). *unit* describes unit of values stored in *point*. Returned value is expressed in *result_unit*.

        This method requires all atomic coordinates to be numerical values, :exc:`~exceptions.TypeError` is raised otherwise.
        """
        ratio = Units.conversion_ratio(unit, 'angstrom')
        resultratio = Units.conversion_ratio('angstrom', result_unit)
        return tuple((i*ratio-j)*resultratio for i,j in zip(point, self))


    def angle(self, point1, point2, point1unit='angstrom', point2unit='angstrom',result_unit='radian'):
        """Calculate an angle between vectors pointing from this atom to *point1* and *point2*.

        *point1* and *point2* should be iterable containers of length 3 (for example: tuple, |Atom|, list, numpy array). Values stored in them are expressed in, respectively, *point1unit* and *point2unit*. Returned value is expressed in *result_unit*.

        This method requires all atomic coordinates to be numerical values, :exc:`~exceptions.TypeError` is raised otherwise.
        """
        num = np.dot(self.vector_to(point1, point1unit), self.vector_to(point2, point2unit))
        den = self.distance_to(point1, point1unit) * self.distance_to(point2, point2unit)
        return Units.convert(math.acos(num/den), 'radian', result_unit)


    def rotate(self, matrix):
        """Rotate this atom according to a rotation *matrix*.

        *matrix* should be a container with 9 numerical values. It can be a list (tuple, numpy array etc.) listing matrix elements row-wise, either flat (``[1,2,3,4,5,6,7,8,9]``) or in two-level fashion (``[[1,2,3],[4,5,6],[7,8,9]]``).

        .. note::

            This method does not check if *matrix* is a proper rotation matrix.
        """
        matrix = np.array(matrix).reshape(3,3)
        self.coords = tuple(np.dot(matrix, np.array(self.coords)))

    def neighbors(self):
        """Return a list of neighbors of this atom within the molecule. The list follows the same order as the ``bonds`` attribute."""
        return [b.other_end(self) for b in self.bonds]
