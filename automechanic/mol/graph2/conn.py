""" connectivity graph module

atm_dct: {atm_key: (atm_sym, atm_hyd_cnt), ...}
bnd_dct: {bnd_key: None, ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
from ._shared import atoms
from ._shared import bonds
from ._shared import atom_keys
from ._shared import bond_keys
from ._shared import atom_symbols
from ._shared import atom_hydrogen_counts
from ._shared import atom_nuclear_charges
from ._shared import atom_total_valences
from ._shared import atom_bonds
from ._shared import atom_neighbor_keys
from ._shared import branch
from ._shared import relabel
from ._shared import subgraph
from ._shared import ring_keys_list
from ._shared import isomorphism
from ._shared import isomorphic
from ._shared import highspin_resonance_graph
from ._shared import backbone_keys
from ._shared import atom_explicit_hydrogen_keys
from ._shared import _from_data
from ._to_inchi import inchi
from ._to_inchi import atom_inchi_numbers


def from_data(atm_keys, bnd_keys, atm_sym_dct, atm_hyd_cnt_dct):
    """ connectivity graph from data
    """
    assert len(atm_keys) == len(atm_sym_dct), len(atm_hyd_cnt_dct)
    return _from_data(
        atm_keys=atm_keys,
        bnd_keys=bnd_keys,
        atm_dcts=[atm_sym_dct, atm_hyd_cnt_dct],
        bnd_dct={}
    )


__all__ = ['atoms', 'bonds', 'atom_keys', 'bond_keys', 'atom_symbols',
           'atom_hydrogen_counts', 'atom_nuclear_charges',
           'atom_total_valences', 'atom_bonds', 'atom_neighbor_keys',
           'atom_explicit_hydrogen_keys',
           'branch', 'relabel', 'subgraph', 'ring_keys_list', 'backbone_keys',
           'isomorphism', 'isomorphic', 'highspin_resonance_graph', 'inchi',
           'atom_inchi_numbers']
