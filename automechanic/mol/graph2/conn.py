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
from ._shared import atom_neighbor_keys
from ._shared import relabel
from ._shared import isomorphism
from ._shared import isomorphic
from ._shared import _from_data


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
           'atom_hydrogen_counts', 'atom_neighbor_keys', 'relabel',
           'isomorphism', 'isomorphic']
