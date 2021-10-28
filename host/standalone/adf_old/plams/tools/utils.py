from __future__ import unicode_literals

import collections
import math
import numpy

from ..core.errors import PTError, UnitsError

__all__ = ['PeriodicTable', 'PT', 'Units']

class PeriodicTable(object):
    """Singleton class for periodic table of elements.

    For each element the following properties are stores: atomic symbol, atomic mass, atomic radius and number of connectors.

    Atomic mass is, strictly speaking, atomic weight, as present in Mathematica's ElementData function.

    Atomic radius and number of connectors are used by :meth:`~scm.plams.basemol.Molecule.guess_bonds`. Note that values or radii are neither atomic radii nor covalent radii. They are someway "emprically optimized" for bond guessing algorithm.

    .. note::

        This class is visible in the main namespace as both ``PeriodicTable`` and ``PT``.
    """
    data = [None] * 113
               #[symbol, mass, radius, connectors]
    data[  0] = ['Xx',   0.00000, 0.00 ,  0]
    data[  1] = [ 'H',   1.00794, 0.30 ,  1]
    data[  2] = ['He',   4.00260, 0.99 ,  0]
    data[  3] = ['Li',   6.94100, 1.52 ,  8]
    data[  4] = ['Be',   9.01218, 1.12 ,  8]
    data[  5] = [ 'B',  10.81100, 0.88 ,  6]
    data[  6] = [ 'C',  12.01070, 0.77 ,  4]
    data[  7] = [ 'N',  14.00670, 0.70 ,  3]
    data[  8] = [ 'O',  15.99940, 0.66 ,  2]
    data[  9] = [ 'F',  18.99840, 0.64 ,  1]
    data[ 10] = ['Ne',  20.17970, 1.60 ,  0]
    data[ 11] = ['Na',  22.98977, 1.86 ,  8]
    data[ 12] = ['Mg',  24.30500, 1.60 ,  8]
    data[ 13] = ['Al',  26.98154, 1.43 ,  8]
    data[ 14] = ['Si',  28.08550, 1.17 ,  8]
    data[ 15] = [ 'P',  30.97376, 1.10 ,  8]
    data[ 16] = [ 'S',  32.06500, 1.04 ,  2]
    data[ 17] = ['Cl',  35.45300, 0.99 ,  1]
    data[ 18] = ['Ar',  39.94800, 1.92 ,  0]
    data[ 19] = [ 'K',  39.09830, 2.31 ,  8]
    data[ 20] = ['Ca',  40.07800, 1.97 ,  8]
    data[ 21] = ['Sc',  44.95591, 1.60 ,  8]
    data[ 22] = ['Ti',  47.86700, 1.46 ,  8]
    data[ 23] = [ 'V',  50.94150, 1.31 ,  8]
    data[ 24] = ['Cr',  51.99610, 1.25 ,  8]
    data[ 25] = ['Mn',  54.93805, 1.29 ,  8]
    data[ 26] = ['Fe',  55.84500, 1.26 ,  8]
    data[ 27] = ['Co',  58.93320, 1.25 ,  8]
    data[ 28] = ['Ni',  58.69340, 1.24 ,  8]
    data[ 29] = ['Cu',  63.54600, 1.28 ,  8]
    data[ 30] = ['Zn',  65.40900, 1.33 ,  8]
    data[ 31] = ['Ga',  69.72300, 1.41 ,  8]
    data[ 32] = ['Ge',  72.64000, 1.22 ,  8]
    data[ 33] = ['As',  74.92160, 1.21 ,  8]
    data[ 34] = ['Se',  78.96000, 1.17 ,  8]
    data[ 35] = ['Br',  79.90400, 1.14 ,  1]
    data[ 36] = ['Kr',  83.79800, 1.97 ,  0]
    data[ 37] = ['Rb',  85.46780, 2.44 ,  8]
    data[ 38] = ['Sr',  87.62000, 2.15 ,  8]
    data[ 39] = [ 'Y',  88.90585, 1.80 ,  8]
    data[ 40] = ['Zr',  91.22400, 1.57 ,  8]
    data[ 41] = ['Nb',  92.90638, 1.41 ,  8]
    data[ 42] = ['Mo',  95.94000, 1.36 ,  8]
    data[ 43] = ['Tc',  98.00000, 1.35 ,  8]
    data[ 44] = ['Ru', 101.07000, 1.33 ,  8]
    data[ 45] = ['Rh', 102.90550, 1.34 ,  8]
    data[ 46] = ['Pd', 106.42000, 1.38 ,  8]
    data[ 47] = ['Ag', 107.86820, 1.44 ,  8]
    data[ 48] = ['Cd', 112.41100, 1.49 ,  8]
    data[ 49] = ['In', 114.81800, 1.66 ,  8]
    data[ 50] = ['Sn', 118.71000, 1.62 ,  8]
    data[ 51] = ['Sb', 121.76000, 1.41 ,  8]
    data[ 52] = ['Te', 127.60000, 1.37 ,  8]
    data[ 53] = [ 'I', 126.90447, 1.33 ,  1]
    data[ 54] = ['Xe', 131.29300, 2.17 ,  0]
    data[ 55] = ['Cs', 132.90545, 2.62 ,  8]
    data[ 56] = ['Ba', 137.32700, 2.17 ,  8]
    data[ 57] = ['La', 138.90550, 1.88 ,  8]
    data[ 58] = ['Ce', 140.11600, 1.818,  8]
    data[ 59] = ['Pr', 140.90765, 1.824,  8]
    data[ 60] = ['Nd', 144.24000, 1.814,  8]
    data[ 61] = ['Pm', 145.00000, 1.834,  8]
    data[ 62] = ['Sm', 150.36000, 1.804,  8]
    data[ 63] = ['Eu', 151.96400, 2.084,  8]
    data[ 64] = ['Gd', 157.25000, 1.804,  8]
    data[ 65] = ['Tb', 158.92534, 1.773,  8]
    data[ 66] = ['Dy', 162.50000, 1.781,  8]
    data[ 67] = ['Ho', 164.93032, 1.762,  8]
    data[ 68] = ['Er', 167.25900, 1.761,  8]
    data[ 69] = ['Tm', 168.93421, 1.759,  8]
    data[ 70] = ['Yb', 173.04000, 1.922,  8]
    data[ 71] = ['Lu', 174.96700, 1.738,  8]
    data[ 72] = ['Hf', 178.49000, 1.57 ,  8]
    data[ 73] = ['Ta', 180.94790, 1.43 ,  8]
    data[ 74] = [ 'W', 183.84000, 1.37 ,  8]
    data[ 75] = ['Re', 186.20700, 1.37 ,  8]
    data[ 76] = ['Os', 190.23000, 1.34 ,  8]
    data[ 77] = ['Ir', 192.21700, 1.35 ,  8]
    data[ 78] = ['Pt', 195.07800, 1.38 ,  8]
    data[ 79] = ['Au', 196.96655, 1.44 ,  8]
    data[ 80] = ['Hg', 200.59000, 1.52 ,  8]
    data[ 81] = ['Tl', 204.38330, 1.71 ,  8]
    data[ 82] = ['Pb', 207.20000, 1.75 ,  8]
    data[ 83] = ['Bi', 208.98038, 1.70 ,  8]
    data[ 84] = ['Po', 209.00000, 1.40 ,  8]
    data[ 85] = ['At', 210.00000, 1.40 ,  1]
    data[ 86] = ['Rn', 222.00000, 2.40 ,  0]
    data[ 87] = ['Fr', 223.00000, 2.70 ,  8]
    data[ 88] = ['Ra', 226.00000, 2.20 ,  8]
    data[ 89] = ['Ac', 227.00000, 2.00 ,  8]
    data[ 90] = ['Th', 232.03810, 1.79 ,  8]
    data[ 91] = ['Pa', 231.03588, 1.63 ,  8]
    data[ 92] = [ 'U', 238.02891, 1.56 ,  8]
    data[ 93] = ['Np', 237.00000, 1.55 ,  8]
    data[ 94] = ['Pu', 244.00000, 1.59 ,  8]
    data[ 95] = ['Am', 243.00000, 1.73 ,  8]
    data[ 96] = ['Cm', 247.00000, 1.74 ,  8]
    data[ 97] = ['Bk', 247.00000, 1.70 ,  8]
    data[ 98] = ['Cf', 251.00000, 1.86 ,  8]
    data[ 99] = ['Es', 252.00000, 1.86 ,  8]
    data[100] = ['Fm', 257.00000, 2.00 ,  8]
    data[101] = ['Md', 258.00000, 2.00 ,  8]
    data[102] = ['No', 259.00000, 2.00 ,  8]
    data[103] = ['Lr', 262.00000, 2.00 ,  8]
    data[104] = ['Rf', 261.00000, 2.00 ,  8]
    data[105] = ['Db', 262.00000, 2.00 ,  8]
    data[106] = ['Sg', 266.00000, 2.00 ,  8]
    data[107] = ['Bh', 264.00000, 2.00 ,  8]
    data[108] = ['Hs', 277.00000, 2.00 ,  8]
    data[109] = ['Mt', 268.00000, 2.00 ,  8]
    data[110] = ['Ds', 281.00000, 2.00 ,  8]
    data[111] = ['Rg', 280.00000, 2.00 ,  8]
    data[112] = ['Cn', 285.00000, 2.00 ,  8]

    symtonum = {d[0]:i for i,d in enumerate(data)}


    def __init__(self):
        raise PTError('Instances of PeriodicTable cannot be created')


    @classmethod
    def get_atomic_number(cls, symbol):
        """Convert atomic symbol to atomic number."""
        try:
            number = cls.symtonum[symbol.capitalize()]
        except KeyError:
            raise PTError('trying to convert incorrect atomic symbol')
        return number


    @classmethod
    def get_symbol(cls, atnum):
        """Convert atomic number to atomic symbol."""
        try:
            symbol = cls.data[atnum][0]
        except IndexError:
            raise PTError('trying to convert incorrect atomic number')
        return symbol


    @classmethod
    def get_mass(cls, arg):
        """Convert atomic symbol or atomic number to atomic mass."""
        return cls._get_property(arg, 1)


    @classmethod
    def get_radius(cls, arg):
        """Convert atomic symbol or atomic number to radius."""
        return cls._get_property(arg, 2)


    @classmethod
    def get_connectors(cls, arg):
        """Convert atomic symbol or atomic number to number of connectors."""
        return cls._get_property(arg, 3)


    @classmethod
    def _get_property(cls, arg, prop):
        """Get property of element described by either symbol or atomic number. Skeleton method for :meth`get_radius`, :meth`get_mass` and  :meth`get_connectors`."""
        if isinstance(arg, str):
            pr = cls.data[cls.get_atomic_number(arg)][prop]
        elif isinstance(arg, int):
            try:
                pr = cls.data[arg][prop]
            except KeyError:
                raise PTError('trying to convert incorrect atomic number')
        return pr



PT = PeriodicTable


#===================================================================================================
#===================================================================================================
#===================================================================================================


class Units(object):
    """Singleton class for units converter.

    All values are based on `2014 CODATA recommended values <http://physics.nist.gov/cuu/Constants>`_.

    The following constants and units are supported:
        *   constants:

            -   ``speed_of_light`` (also ``c``)
            -   ``elementary_charge`` (also ``e`` and ``electron_charge``)
            -   ``avogadro_constant`` (also ``NA``)
            -   ``bohr_radius``

        *   distance:

            -   ``Angstrom``, ``angstrom``, ``A``
            -   ``bohr``, ``a0``, ``au``
            -   ``nm``
            -   ``pm``

        *   angle:

            -    ``degree``, ``deg``,
            -    ``radian``, ``rad``,
            -    ``grad``
            -    ``circle``

        *   energy:

            -   ``au``, ``hartree``, ``Hartree``
            -   ``ev``, ``eV``
            -   ``kcal/mol``
            -   ``kJ/mol``
            -   ``cm^-1``

        *   dipole moment:

            -   ``au``
            -   ``Cm``
            -   ``D``, ``Debye``, ``debye``

    Example::

        >>> print(Units.constants['speed_of_light'])
        299792458
        >>> print(Units.constants['e'])
        1.6021766208e-19
        >>> print(Units.convert(123, 'angstrom', 'bohr'))
        232.436313431
        >>> print(Units.convert(23.32, 'kJ/mol', 'kcal/mol'))
        5.57361376673
        >>> print(Units.conversion_ratio('kcal/mol', 'kJ/mol'))
        4.184


    """

    constants = {}
    constants['bohr_radius'] = 0.52917721067   #http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
    constants['avogadro_constant'] = 6.022140857e23   #http://physics.nist.gov/cgi-bin/cuu/Value?na
    constants['speed_of_light'] = 299792458   #http://physics.nist.gov/cgi-bin/cuu/Value?c
    constants['electron_charge'] = 1.6021766208e-19   #http://physics.nist.gov/cgi-bin/cuu/Value?e

    constants['NA'] = constants['avogadro_constant']
    constants['c'] = constants['speed_of_light']
    constants['e'] = constants['electron_charge']
    constants['elementary_charge'] = constants['electron_charge']

    dicts = []

    distance = {}
    distance['A'] = 1.0
    distance['angstrom'] = distance['A']
    distance['Angstrom'] = distance['A']
    distance['nm'] = distance['A'] / 10.0
    distance['pm'] = distance['A'] * 100.0
    distance['bohr'] = 1.0 / constants['bohr_radius']
    distance['a0'] = distance['bohr']
    distance['au'] = distance['bohr']
    dicts.append(distance)

    energy = {}
    energy['au'] = 1.0
    energy['hartree'] = energy['au']
    energy['Hartree'] = energy['au']
    energy['eV'] = 27.21138602   #http://physics.nist.gov/cgi-bin/cuu/Value?hrev
    energy['ev'] = energy['eV']
    energy['kJ/mol'] = 4.359744650e-21 * constants['NA']  #http://physics.nist.gov/cgi-bin/cuu/Value?hrj
    energy['kcal/mol'] = energy['kJ/mol'] / 4.184
    energy['cm^-1'] = 219474.6313702   #http://physics.nist.gov/cgi-bin/cuu/Value?hrminv
    dicts.append(energy)

    angle = {}
    angle['degree'] = 1.0
    angle['deg'] = angle['degree']
    angle['radian'] = math.pi / 180.0
    angle['rad'] = angle['radian']
    angle['grad'] = 100.0 / 90.0
    angle['circle'] = 1.0 / 360.0
    dicts.append(angle)

    dipole = {}
    dipole['au'] = 1.0
    dipole['Cm'] = constants['e'] * constants['bohr_radius'] * 1e-10
    dipole['Debye'] = dipole['Cm'] * constants['c']* 1e21
    dipole['debye'] = dipole['Debye']
    dipole['D'] = dipole['Debye']
    dicts.append(dipole)

    def __init__(self):
        raise UnitsError('Instances of Units cannot be created')


    @classmethod
    def conversion_ratio(cls, inp, out):
        """Return conversion ratio from unit *inp* to *out*."""
        for d in cls.dicts:
            if inp in d.keys() and out in d.keys():
                return d[out]/d[inp]
        raise UnitsError('Invalid conversion_ratio call: unsupported units')


    @classmethod
    def convert(cls, value, inp, out):
        """Convert *value* from unit *inp* to *out*. *value* has to be either a single number or a container (list, tuple etc.) containing only numbers. In this case container of the same type and length is returned."""
        r = cls.conversion_ratio(inp,out)
        if isinstance(value, collections.Iterable):
            t = type(value)
            if t == numpy.ndarray:
                t = numpy.array
            v = [i*r for i in value]
            return t(v)
        return value * r
