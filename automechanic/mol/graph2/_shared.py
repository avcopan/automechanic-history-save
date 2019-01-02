""" shared molecular graph functions
"""
from ._inetworkx import from_graph as _nxg_from_graph
from ._inetworkx import isomorphism as _nxg_isomorphism

SYMBOL_POS = 0
HCOUNT_POS = 1
PARITY_POS = 2


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


def _atom_values_at_position(xgr, pos):
    atm_dct = atoms(xgr)
    atm_keys = atm_dct.keys()
    atm_vals = atm_dct.values()
    pos_atm_vals = list(zip(*atm_vals))[pos]
    return dict(zip(atm_keys, pos_atm_vals))


def atom_symbols(xgr):
    """ atom symbols, as a dictionary
    """
    return _atom_values_at_position(xgr, SYMBOL_POS)


def atom_hydrogen_counts(xgr):
    """ atom symbols, as a dictionary
    """
    return _atom_values_at_position(xgr, HCOUNT_POS)


def _atom_stereo_parities(xgr):
    return _atom_values_at_position(xgr, PARITY_POS)


def atom_neighbor_keys(xgr):
    """ keys of neighboring atoms, by atom
    """
    atm_keys = atom_keys(xgr)
    bnd_keys = bond_keys(xgr)
    return {atm_key: tuple(sorted(next(iter(bnd_key - {atm_key}))
                                  for bnd_key in bnd_keys
                                  if atm_key in bnd_key))
            for atm_key in atm_keys}


def relabel(xgr, atm_key_dct):
    """ relabel the graph with new atom keys
    """
    atm_dct = atoms(xgr)
    bnd_dct = bonds(xgr)
    atm_keys = atm_dct.keys()
    atm_vals = atm_dct.values()
    bnd_keys = bnd_dct.keys()
    bnd_vals = bnd_dct.values()

    assert set(atm_key_dct.keys()) == set(atom_keys(xgr))

    atm_keys = list(map(atm_key_dct.__getitem__, atm_keys))
    bnd_keys = list(map(
        frozenset,
        (map(atm_key_dct.__getitem__, bnd_key) for bnd_key in bnd_keys)))

    atm_dct = dict(zip(atm_keys, atm_vals))
    bnd_dct = dict(zip(bnd_keys, bnd_vals))
    xgr = (atm_dct, bnd_dct)
    return xgr


def isomorphic(xgr1, xgr2):
    """ are these molecular graphs isomorphic?
    """
    return isomorphism(xgr1, xgr2) is not None


def isomorphism(xgr1, xgr2):
    """ graph isomorphism (relabeling of `xgr1` to produce `xgr2`)
    """
    nxg1 = _nxg_from_graph(xgr1)
    nxg2 = _nxg_from_graph(xgr2)
    iso_dct = _nxg_isomorphism(nxg1, nxg2)
    return iso_dct


def _from_data(atm_keys, bnd_keys, atm_dcts, bnd_dct):
    atm_dcts = [_fill_nones(atm_dct, atm_keys) for atm_dct in atm_dcts]
    bnd_dct = _fill_nones(bnd_dct, bnd_keys)

    assert all(set(atm_dct.keys()) == set(atm_keys) for atm_dct in atm_dcts)
    assert set(bnd_dct.keys()) == set(bnd_keys)

    atm_dct = dict(zip(
        atm_keys,
        zip(*(map(atm_dct.__getitem__, atm_keys) for atm_dct in atm_dcts))
    ))

    xgr = (atm_dct, bnd_dct)
    return xgr


def _fill_nones(dct, keys):
    return dict((key, dct[key]) if key in dct else (key, None) for key in keys)
