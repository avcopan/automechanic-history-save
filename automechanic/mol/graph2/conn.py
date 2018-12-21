""" connectivity graph module

atms: frozenset({(atm_key, (atm_sym, atm_hyd_cnt)), ...})
bnds: frozenset({(bnd_key, None), ...})
bnd_key := frozenset({atm1_key, atm2_key})
"""
from ._shared import atom_keys
from ._shared import atom_symbols
from ._shared import atom_hydrogen_counts
from ._shared import _from_data
from ._shared import bond_keys


def from_data(atm_keys, atm_syms, atm_hyd_cnts, bnd_keys):
    """ connectivity graph from data
    """
    assert len(atm_keys) == len(atm_syms) == len(atm_hyd_cnts)

    atm_vals = tuple(zip(atm_syms, atm_hyd_cnts))
    bnd_vals = (None,) * len(bnd_keys)
    return _from_data(atm_keys=atm_keys, atm_vals=atm_vals,
                      bnd_keys=bnd_keys, bnd_vals=bnd_vals)


__all__ = ['atom_keys', 'atom_symbols', 'atom_hydrogen_counts', 'bond_keys']
