""" molecular graph modules

xgr = (atm_dct, bnd_dct)
atm_dct: {atm_key: (atm_val1, atm_val2, ...), ...}
bnd_dct: {bnd_key: (bnd_val1, bnd_val2, ...), ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
from . import res
from . import stereo
from ._base import from_data
from ._base import atoms
from ._base import bonds
from ._base import atom_keys
from ._base import bond_keys
from ._base import atom_symbols
from ._base import atom_implicit_hydrogen_valences
from ._base import set_atom_implicit_hydrogen_valences

__all__ = ['res', 'stereo',
           'from_data', 'atoms', 'bonds', 'atom_keys', 'bond_keys',
           'atom_symbols', 'atom_implicit_hydrogen_valences',
           'set_atom_implicit_hydrogen_valences']
