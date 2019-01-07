""" base graph library; depending only on connectivity

xgr = (atm_dct, bnd_dct)
atm_dct: {atm_key: (atm_sym, atm_imp_hyd_vlc, ...), ...}
bnd_dct: {bnd_key: (1 or bnd_ord, ...), ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
from ._dict import by_key as _by_key
from ._dict import values_by_key as _values_by_key
from ._dict import transform_values as _transform_values
from ._dict import transform_values_with_key as _transform_values_with_key
from ._dict import filter_by_key as _filter_by_key
from ._dict import filter_by_value as _filter_by_value
from ._tdict import by_key_by_position as _by_key_by_position
from ._tdict import set_by_key_by_position as _set_by_key_by_position
from ..atom import nuclear_charge as _atom_nuclear_charge
from ..atom import valence as _atom_valence

ATM_SYM_POS = 0
ATM_IMP_HYD_VLC_POS = 1

# not used here, but for reference:
# ATM_PAR_POS = 2
# BND_ORD_POS = 0
# BND_PAR_POS = 1


def from_data(atm_keys, bnd_keys, atm_sym_dct, atm_imp_hyd_vlc_dct=None,
              atm_par_dct=None, bnd_ord_dct=None, bnd_par_dct=None):
    """ molecular graph (any type) from data
    """
    atm_val_dcts = []
    bnd_val_dcts = []

    # check each contribution and append it to the list
    assert set(atm_sym_dct.keys()) == set(atm_keys)
    atm_val_dcts.append(atm_sym_dct)

    atm_imp_hyd_vlc_dct = ({} if atm_imp_hyd_vlc_dct is None else
                           atm_imp_hyd_vlc_dct)
    assert set(atm_imp_hyd_vlc_dct.keys()) <= set(atm_keys)
    atm_imp_hyd_vlc_dct = _by_key(atm_imp_hyd_vlc_dct, atm_keys, fill_val=0)
    atm_val_dcts.append(atm_imp_hyd_vlc_dct)

    bnd_ord_dct = {} if bnd_ord_dct is None else bnd_ord_dct
    assert set(bnd_ord_dct.keys()) <= set(bnd_keys)
    bnd_ord_dct = _by_key(bnd_ord_dct, bnd_keys, fill_val=1)
    bnd_val_dcts.append(bnd_ord_dct)

    if atm_par_dct is not None or bnd_par_dct is not None:
        atm_par_dct = {} if atm_par_dct is None else atm_par_dct
        assert set(atm_par_dct.keys()) <= set(atm_keys)
        atm_val_dcts.append(atm_par_dct)

        bnd_par_dct = {} if bnd_par_dct is None else bnd_par_dct
        assert set(bnd_par_dct.keys()) <= set(bnd_keys)
        bnd_val_dcts.append(bnd_par_dct)

    return _from_data(atm_keys, bnd_keys, atm_val_dcts, bnd_val_dcts)


def atoms(xgr):
    """ atoms, as a dictionary
    """
    atm_dct, _ = xgr
    return atm_dct


def bonds(xgr):
    """ bonds, as a dictionary
    """
    _, bnd_dct = xgr
    return bnd_dct


def atom_keys(xgr):
    """ sorted atom keys
    """
    return tuple(sorted(atoms(xgr).keys()))


def bond_keys(xgr):
    """ sorted bond keys
    """
    return tuple(sorted(bonds(xgr).keys(), key=sorted))


def atom_symbols(xgr):
    """ atom symbols, as a dictionary
    """
    return _by_key_by_position(atoms(xgr), atom_keys(xgr), ATM_SYM_POS)


def atom_implicit_hydrogen_valences(xgr):
    """ atom implicit hydrogen valences, as a dictionary
    """
    return _filter_by_value(
        _by_key_by_position(atoms(xgr), atom_keys(xgr), ATM_IMP_HYD_VLC_POS))


def set_atom_implicit_hydrogen_valences(xgr, atm_imp_hyd_vlc_dct):
    """ set atom implicit hydrogen valences
    """
    atm_dct = _set_by_key_by_position(atoms(xgr), atm_imp_hyd_vlc_dct,
                                      ATM_IMP_HYD_VLC_POS)
    bnd_dct = bonds(xgr)
    xgr = (atm_dct, bnd_dct)
    return xgr


def atom_nuclear_charges(xgr):
    """ nuclear charges, by atom (connectivity-independent)
    """
    atm_sym_dct = atom_symbols(xgr)
    atm_nuc_chg_dct = _transform_values(atm_sym_dct, func=_atom_nuclear_charge)
    return atm_nuc_chg_dct


def atom_total_valences(xgr):
    """ total valence electron counts, by atom (connectivity-independent)
    """
    atm_sym_dct = atom_symbols(xgr)
    atm_tot_vlc_dct = _transform_values(atm_sym_dct, func=_atom_valence)
    return atm_tot_vlc_dct


def atom_bonds(xgr):
    """ bonds, by atom
    """
    def _is_my_bond(atm_key):
        def __is_my_bond(bnd_key):
            return atm_key in bnd_key
        return __is_my_bond

    bnd_dct = bonds(xgr)
    atm_bnds_dct = {atm_key: _filter_by_key(bnd_dct, func=_is_my_bond(atm_key))
                    for atm_key in atom_keys(xgr)}
    return _filter_by_value(atm_bnds_dct)


def atom_neighbor_keys(xgr):
    """ keys of neighboring atoms, by atom
    """
    def _neighbor_keys(atm_key, atm_bnd_dct):
        bnd_keys = atm_bnd_dct.keys()
        ngb_keys = sorted(next(iter(bnd_key - {atm_key}))
                          for bnd_key in bnd_keys)
        return tuple(ngb_keys)

    return _transform_values_with_key(atom_bonds(xgr), func=_neighbor_keys)


# change this to `atom_explicit_hydrogen_valences`
# with separate `backbone_keys` and `explicit_hydrogen_keys` functions
# def atom_explicit_hydrogen_keys(xgr):
#     """ explicit hydrogen keys, by atom
#     """
#     atm_sym_dct = atom_symbols(xgr)
#
#     def _hydrogen_keys(atm_keys):
#         return tuple(atm_key for atm_key in atm_keys
#                      if atm_sym_dct[atm_key] == 'H')
#
#     return _filter_by_value(
#         _transform_values(atom_neighbor_keys(xgr), func=_hydrogen_keys))


def _from_data(atm_keys, bnd_keys, atm_val_dcts, bnd_val_dcts):
    assert all(set(atm_val_dct.keys()) <= set(atm_keys)
               for atm_val_dct in atm_val_dcts)
    assert all(set(bnd_val_dct.keys()) <= set(bnd_keys)
               for bnd_val_dct in bnd_val_dcts)

    atm_val_lsts = [_values_by_key(atm_val_dct, atm_keys)
                    for atm_val_dct in atm_val_dcts]
    bnd_val_lsts = [_values_by_key(bnd_val_dct, bnd_keys)
                    for bnd_val_dct in bnd_val_dcts]

    atm_dct = dict(zip(atm_keys, zip(*atm_val_lsts)))
    bnd_dct = dict(zip(bnd_keys, zip(*bnd_val_lsts)))

    xgr = (atm_dct, bnd_dct)
    return xgr