""" stereo graph module

atms: frozenset({(atm_key, (atm_sym, atm_hyd_cnt, atm_par)), ...})
bnds: frozenset({(bnd_key, bnd_par), ...})
bnd_key := frozenset({atm1_key, atm2_key})
"""
from ._shared import atom_keys
from ._shared import atom_symbols
from ._shared import atom_hydrogen_counts
from ._shared import atom_stereo_parities
from ._shared import bond_keys
from ._shared import _from_data
from ._shared import _bond_values


def bond_stereo_parities(rgr):
    """ bond stereo parities, sorted by bond key
    """
    return _bond_values(rgr)


def from_data(atm_keys, atm_syms, atm_hyd_cnts, atm_pars, bnd_keys, bnd_pars):
    """ resonance graph from data
    """
    assert len(atm_keys) == len(atm_syms) == len(atm_hyd_cnts) == len(atm_pars)
    assert len(bnd_keys) == len(bnd_pars)

    atm_vals = tuple(zip(atm_syms, atm_hyd_cnts, atm_pars))
    bnd_vals = bnd_pars
    return _from_data(atm_keys=atm_keys, atm_vals=atm_vals,
                      bnd_keys=bnd_keys, bnd_vals=bnd_vals)


__all__ = ['atom_keys', 'atom_symbols', 'atom_hydrogen_counts',
           'atom_stereo_parities', 'bond_keys']
