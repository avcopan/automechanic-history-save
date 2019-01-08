""" molecular graph modules

xgr = (atm_dct, bnd_dct)
atm_dct: {atm_key: (atm_val1, atm_val2, ...), ...}
bnd_dct: {bnd_key: (bnd_val1, bnd_val2, ...), ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
from . import res
from . import stereo
from ._base import add_atom
from ._base import from_data
from ._base import atoms
from ._base import bonds
from ._base import atom_keys
from ._base import bond_keys
from ._base import atom_symbols
from ._base import atom_implicit_hydrogen_valences
from ._base import set_atom_implicit_hydrogen_valences
from ._base import atom_nuclear_charges
from ._base import atom_total_valences
from ._base import atom_bonds
from ._base import atom_neighbor_keys
from ._base import explicit_hydrogen_keys
from ._base import backbone_keys
from ._base import atom_explicit_hydrogen_keys
from ._base import implicit
from ._base import delete_atoms
from ._base import subgraph

__all__ = ['res', 'stereo',
           'add_atom',
           'from_data', 'atoms', 'bonds', 'atom_keys', 'bond_keys',
           'atom_symbols', 'atom_implicit_hydrogen_valences',
           'set_atom_implicit_hydrogen_valences', 'atom_nuclear_charges',
           'atom_total_valences', 'atom_bonds', 'atom_neighbor_keys',
           'explicit_hydrogen_keys', 'backbone_keys',
           'atom_explicit_hydrogen_keys', 'delete_atoms', 'implicit',
           'subgraph']
