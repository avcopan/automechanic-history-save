""" resonance graph module

atm_dct: {atm_key: (atm_sym, atm_hyd_cnt), ...}
bnd_dct: {bnd_key: bnd_cnt, ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
import numpy
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
from ._shared import isomorphism
from ._shared import isomorphic
from ._shared import highspin_resonance_graph
from ._shared import _from_data
from ._dict import values_by_key as _values_by_key


def from_data(atm_keys, bnd_keys, atm_sym_dct, atm_hyd_cnt_dct, bnd_ord_dct):
    """ resonance graph from data
    """
    assert len(atm_keys) == len(atm_sym_dct), len(atm_hyd_cnt_dct)
    assert len(bnd_keys) == len(bnd_ord_dct)
    return _from_data(
        atm_keys=atm_keys,
        bnd_keys=bnd_keys,
        atm_dcts=[atm_sym_dct, atm_hyd_cnt_dct],
        bnd_dct=bnd_ord_dct
    )


def bond_orders(rgr):
    """ bond orders (alias of bonds)
    """
    return bonds(rgr)


def atom_bond_valences(rgr):
    """ bond valences, by atom
    """
    atm_keys = atom_keys(rgr)
    atm_hyd_cnts = _values_by_key(atom_hydrogen_counts(rgr), atm_keys)
    atm_bnd_dcts = _values_by_key(atom_bonds(rgr), atm_keys)
    atm_bnd_cnts = [sum(atm_bnd_dct.values()) for atm_bnd_dct in atm_bnd_dcts]
    atm_bnd_vlcs = numpy.add(atm_hyd_cnts, atm_bnd_cnts)
    return dict(zip(atm_keys, atm_bnd_vlcs))


def atom_radical_valences(rgr):
    """ radical valences, by atom
    """
    atm_keys = atom_keys(rgr)
    atm_bnd_vlcs = _values_by_key(atom_bond_valences(rgr), atm_keys)
    atm_tot_vlcs = _values_by_key(atom_total_valences(rgr), atm_keys)
    atm_rad_vlcs = numpy.subtract(atm_tot_vlcs, atm_bnd_vlcs)
    return {atm_key: atm_rad_vlc
            for atm_key, atm_rad_vlc in zip(atm_keys, atm_rad_vlcs)
            if atm_rad_vlc}


__all__ = ['atoms', 'bonds', 'atom_keys', 'bond_keys', 'atom_symbols',
           'atom_hydrogen_counts', 'atom_nuclear_charges',
           'atom_total_valences', 'atom_bonds', 'atom_neighbor_keys',
           'branch', 'relabel', 'subgraph',
           'isomorphism', 'isomorphic', 'highspin_resonance_graph']
