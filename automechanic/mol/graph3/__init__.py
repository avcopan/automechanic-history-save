""" molecular graph modules

xgr = (atm_dct, bnd_dct)
atm_dct: {atm_key: (atm_val1, atm_val2, ...), ...}
bnd_dct: {bnd_key: (bnd_val1, bnd_val2, ...), ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
# constructors
from ._base import empty_graph
from ._base import from_data
from ._base import add_atoms
from ._base import add_bonds
# value getters
from ._base import atoms
from ._base import bonds
from ._base import atom_keys
from ._base import bond_keys
from ._base import atom_symbols
from ._base import atom_implicit_hydrogen_valences
from ._base import atom_stereo_keys
from ._base import atom_stereo_parities
from ._base import bond_orders
from ._base import bond_stereo_keys
from ._base import bond_stereo_parities
# value setters
from ._base import set_atom_implicit_hydrogen_valences
from ._base import set_atom_stereo_parities
from ._base import set_bond_orders
from ._base import set_bond_stereo_parities
# derived values
from ._base import atom_nuclear_charges
from ._base import atom_bonds
from ._base import atom_neighbor_keys
from ._base import explicit_hydrogen_keys
from ._base import backbone_keys
from ._base import atom_explicit_hydrogen_keys
from ._base import ring_keys_list
from ._base import atom_total_valences
from ._base import atom_bond_valences
from ._base import atom_radical_valences
from ._base import is_chiral
from ._inchi.conn import atom_inchi_numbers
from ._inchi.conn import inchi
from ._inchi.stereo import inchi as stereo_inchi
# transformations
from ._base import implicit
from ._base import explicit
from ._base import explicit_stereo_sites
from ._base import delete_atoms
from ._base import add_explicit_hydrogens
from ._base import subgraph
from ._base import relabel
from ._base import reflection
# comparisons
from ._base import backbone_isomorphic
from ._base import backbone_isomorphism


__all__ = [
    # constructors
    'empty_graph', 'from_data', 'add_atoms', 'add_bonds',
    # value getters
    'atoms', 'bonds', 'atom_keys', 'bond_keys', 'atom_symbols',
    'atom_implicit_hydrogen_valences', 'atom_stereo_keys',
    'atom_stereo_parities', 'bond_orders', 'bond_stereo_keys',
    'bond_stereo_parities',
    # value setters
    'set_atom_implicit_hydrogen_valences', 'set_atom_stereo_parities',
    'set_bond_orders', 'set_bond_stereo_parities',
    # derived values
    'atom_nuclear_charges', 'atom_bonds', 'atom_neighbor_keys',
    'explicit_hydrogen_keys', 'backbone_keys', 'atom_explicit_hydrogen_keys',
    'ring_keys_list', 'atom_total_valences', 'atom_bond_valences',
    'atom_radical_valences', 'is_chiral', 'atom_inchi_numbers', 'inchi',
    'stereo_inchi',
    # transformations
    'implicit', 'explicit', 'explicit_stereo_sites', 'delete_atoms',
    'add_explicit_hydrogens', 'subgraph', 'relabel', 'reflection',
    # comparisons
    'backbone_isomorphic', 'backbone_isomorphism',
]
