""" resonance graph library

rgr = (atm_dct, bnd_dct)
atm_dct: {atm_key: (atm_sym, atm_imp_hyd_vlc, ...), ...}
bnd_dct: {bnd_key: (bnd_ord, ...), ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
from ._base import atoms as _atoms
from ._base import bonds as _bonds
from ._base import atom_keys as _atom_keys
from ._base import bond_keys as _bond_keys
from ._base import atom_symbols as _atom_symbols
from ._base import (atom_implicit_hydrogen_valences as
                    _atom_implicit_hydrogen_valences)
from ._base import from_data as _from_data
from ._tdict import by_key_by_position as _by_key_by_position
from ._tdict import set_by_key_by_position as _set_by_key_by_position

BND_ORD_POS = 0


def no_pi(xgr):
    """ resonance graph with no pi bonds (= connectivity graph)

    (can be called on an molecular graph type)
    """
    return _from_data(
        _atom_keys(xgr), _bond_keys(xgr), _atom_symbols(xgr),
        atm_imp_hyd_vlc_dct=_atom_implicit_hydrogen_valences(xgr),
        bnd_ord_dct={}
    )


def bond_orders(rgr):
    """ bond orders, as a dictionary
    """
    return _by_key_by_position(_bonds(rgr), _bond_keys(rgr), BND_ORD_POS)


def set_bond_orders(rgr, bnd_ord_dct):
    """ set bond orders
    """
    atm_dct = _atoms(rgr)
    bnd_dct = _set_by_key_by_position(_bonds(rgr), bnd_ord_dct, BND_ORD_POS)
    rgr = (atm_dct, bnd_dct)
    return rgr
