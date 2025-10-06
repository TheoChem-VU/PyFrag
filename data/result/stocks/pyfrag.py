
#-----------------------------------------------------------------------------
# Boilerplate
#-----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function, unicode_literals
import logging
import os
log = logging.getLogger(__name__)

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

# Standard library imports
import csv

# External imports

# Bokeh imports

#-----------------------------------------------------------------------------
# Globals and constants
#-----------------------------------------------------------------------------

__all__ = (
    'PYFRAG',
)

#-----------------------------------------------------------------------------
# General API
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Dev API
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Private API
#-----------------------------------------------------------------------------


def _read_data(name):
    '''
    More data options can be added for instance:
    data = {
    'overlap': [],
    ..............
    }
    and
    for row in reader:
        overlap, irc, ................
        data['overlap'].append(overlap)
    ..............
    '''
    currentPath = os.getcwd()
    filename  = os.path.join(currentPath, str('stocks/PYFRAG.csv'))

    data = {
        'irc' : [],
        'elstat' : [],
        'energytotal' : [],
        'inter' : [],
        'oi' : [],
        'pauli' : [],
        'straintotal' : [],
        'bondlength': [],
        'frag1' : [],
        'frag2' : [],
    }
    with open(filename, 'r') as f:
        next(f)
        reader = csv.reader(f, delimiter=str(','))
        for row in reader:
            irc, elstat, energytotal, inter, oi, pauli, straintotal, bondlength, frag1, frag2 = row
            data['irc'].append(irc)
            data['elstat'].append(float(elstat))
            data['energytotal'].append(float(energytotal))
            data['inter'].append(float(inter))
            data['oi'].append(float(oi))
            data['pauli'].append(float(pauli))
            data['straintotal'].append(float(straintotal))
            data['bondlength'].append(float(bondlength))
            data['frag1'].append(float(frag1))
            data['frag2'].append(float(frag2))
    return data

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

PYFRAG = _read_data('PYFRAG')

