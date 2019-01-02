""" molecular graph modules

atm_dct: {atm_key: (atm_val1, atm_val2), ...}
bnd_dct: {bnd_key: bnd_val, ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
from . import conn
from . import res
from . import stereo

__all__ = ['conn', 'res', 'stereo']
