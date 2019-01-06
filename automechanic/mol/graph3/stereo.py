""" stereo graph library

sgr = (atm_dct, bnd_dct)
atm_dct: {atm_key: (atm_sym, atm_imp_hyd_vlc, ...), ...}
bnd_dct: {bnd_key: (1 or bnd_ord, ...), ...}
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
from .res import bond_orders as _bond_orders

ATM_PAR_POS = 2
BND_PAR_POS = 1


def no_centers(xgr):
    """ stereo graph with no stereo centers

    (can be called on an molecular graph type)
    """
    return _from_data(
        _atom_keys(xgr), _bond_keys(xgr), _atom_symbols(xgr),
        atm_imp_hyd_vlc_dct=_atom_implicit_hydrogen_valences(xgr),
        atm_par_dct={}, bnd_ord_dct=_bond_orders(xgr), bnd_par_dct={}
    )


def atom_parities(sgr):
    """ atom parities, as a dictionary
    """
    return _by_key_by_position(_atoms(sgr), _atom_keys(sgr), ATM_PAR_POS)


def bond_parities(sgr):
    """ bond parities, as a dictionary
    """
    return _by_key_by_position(_bonds(sgr), _bond_keys(sgr), BND_PAR_POS)


def set_atom_parities(sgr, atm_par_dct):
    """ set atom parities
    """
    atm_dct = _set_by_key_by_position(_atoms(sgr), atm_par_dct, ATM_PAR_POS)
    bnd_dct = _bonds(sgr)
    sgr = (atm_dct, bnd_dct)
    return sgr


def set_bond_parities(sgr, bnd_par_dct):
    """ set bond parities
    """
    atm_dct = _atoms(sgr)
    bnd_dct = _set_by_key_by_position(_bonds(sgr), bnd_par_dct, BND_PAR_POS)
    sgr = (atm_dct, bnd_dct)
    return sgr
