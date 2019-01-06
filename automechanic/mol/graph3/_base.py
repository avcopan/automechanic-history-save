""" base graph library; depending only on connectivity

xgr = (atm_dct, bnd_dct)
atm_dct: {atm_key: (atm_sym, atm_imp_hyd_vlc, ...), ...}
bnd_dct: {bnd_key: (1 or bnd_ord, ...), ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
from ._dict import by_key as _by_key
from ._dict import values_by_key as _values_by_key
from ._tdict import by_key_by_position as _by_key_by_position
from ._tdict import set_by_key_by_position as _set_by_key_by_position

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
    return _by_key_by_position(atoms(xgr), atom_keys(xgr), ATM_IMP_HYD_VLC_POS)


def set_atom_implicit_hydrogen_valences(xgr, atm_imp_hyd_vlc_dct):
    """ set atom implicit hydrogen valences
    """
    atm_dct = _set_by_key_by_position(atoms(xgr), atm_imp_hyd_vlc_dct,
                                      ATM_IMP_HYD_VLC_POS)
    bnd_dct = bonds(xgr)
    xgr = (atm_dct, bnd_dct)
    return xgr


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
