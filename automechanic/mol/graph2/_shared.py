""" shared molecular graph functions
"""
from itertools import chain as _chain
import numpy
from ._inetworkx import from_graph as _nxg_from_graph
from ._inetworkx import ring_keys_list as _nxg_ring_keys_list
from ._inetworkx import isomorphism as _nxg_isomorphism
from ._dict import by_key as _by_key
from ._dict import values_by_key as _values_by_key
from ._dict import transform_values as _transform_values
from ..atom import valence as _atom_valence
from ..atom import nuclear_charge as _atom_nuclear_charge

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


def _change_atom_values_at_position(xgr, pos, val_dct):
    atm_dct = atoms(xgr)
    bnd_dct = bonds(xgr)
    atm_dct = _transform_values(atm_dct, func=list)
    assert set(val_dct.keys()) <= set(atm_dct.keys())
    for key, val in val_dct.items():
        atm_dct[key][pos] = val
    atm_dct = _transform_values(atm_dct, func=tuple)
    xgr = (atm_dct, bnd_dct)
    return xgr


def atom_symbols(xgr):
    """ atom symbols, as a dictionary
    """
    return _atom_values_at_position(xgr, SYMBOL_POS)


def atom_hydrogen_counts(xgr):
    """ atom hydrogen counts, as a dictionary
    """
    return _atom_values_at_position(xgr, HCOUNT_POS)


def change_atom_hydrogen_counts(xgr, atm_hyd_cnt_dct):
    """ set the atom hydrogen counts of this molecular graph
    """
    return _change_atom_values_at_position(xgr, HCOUNT_POS, atm_hyd_cnt_dct)


def _atom_stereo_parities(xgr):
    return _atom_values_at_position(xgr, PARITY_POS)


def atom_nuclear_charges(xgr):
    """ nuclear charges, by atom (connectivity-independent)
    """
    atm_keys = atom_keys(xgr)
    atm_syms = _values_by_key(atom_symbols(xgr), atm_keys)
    atm_nuc_chgs = tuple(map(_atom_nuclear_charge, atm_syms))
    return dict(zip(atm_keys, atm_nuc_chgs))


def atom_total_valences(xgr):
    """ total valence electron counts, by atom (connectivity-independent)
    """
    atm_keys = atom_keys(xgr)
    atm_syms = _values_by_key(atom_symbols(xgr), atm_keys)
    atm_vlcs = tuple(map(_atom_valence, atm_syms))
    return dict(zip(atm_keys, atm_vlcs))


def atom_bonds(xgr):
    """ bonds, by atom
    """
    atm_keys = atom_keys(xgr)
    bnds = bonds(xgr)
    return {atm_key: {bnd_key: bnd_val for bnd_key, bnd_val in bnds.items()
                      if atm_key in bnd_key}
            for atm_key in atm_keys}


def atom_neighbor_keys(xgr):
    """ keys of neighboring atoms, by atom
    """
    atm_keys = atom_keys(xgr)
    bnd_keys = bond_keys(xgr)
    return {atm_key: tuple(sorted(next(iter(bnd_key - {atm_key}))
                                  for bnd_key in bnd_keys
                                  if atm_key in bnd_key))
            for atm_key in atm_keys}


def branch_keys(xgr, atm_key, atm_ngb_key, excluded_atm_keys=None):
    """ keys for the branch extending from an atom towards one of its neighbors

    exclude_atm_keys can be used to cut off cycles
    """
    assert frozenset({atm_key, atm_ngb_key}) in bond_keys(xgr)
    excluded_atm_keys = (set([]) if excluded_atm_keys is None else
                         set(excluded_atm_keys))

    atm_ngb_keys_dct = atom_neighbor_keys(xgr)
    bnch_atm_keys = {atm_ngb_key}
    excl_atm_keys = {atm_key} | excluded_atm_keys
    seen_atm_keys = set([])

    def _branch_keys_recursive(bnch_atm_keys):
        new_atm_keys = bnch_atm_keys - seen_atm_keys
        if new_atm_keys:
            new_atm_ngb_keys = set(_chain(
                *_values_by_key(atm_ngb_keys_dct, new_atm_keys)))
            new_atm_ngb_keys -= excl_atm_keys
            bnch_atm_keys.update(new_atm_ngb_keys)
            seen_atm_keys.update(new_atm_keys)
            _branch_keys_recursive(bnch_atm_keys)
        return bnch_atm_keys

    bnch_atm_keys = _branch_keys_recursive(bnch_atm_keys)
    bnch_atm_keys.add(atm_key)

    return tuple(sorted(bnch_atm_keys))


def branch(xgr, atm_key, atm_ngb_key, excluded_atm_keys=None):
    """ keys for the branch extending from an atom towards one of its neighbors

    exclude_atm_keys can be used to cut off cycles
    """
    bnch_atm_keys = branch_keys(xgr, atm_key, atm_ngb_key, excluded_atm_keys)
    return subgraph(xgr, bnch_atm_keys)


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

    atm_keys = _values_by_key(atm_key_dct, atm_keys)
    bnd_keys = list(map(
        frozenset,
        (_values_by_key(atm_key_dct, bnd_key) for bnd_key in bnd_keys)))

    atm_dct = dict(zip(atm_keys, atm_vals))
    bnd_dct = dict(zip(bnd_keys, bnd_vals))
    xgr = (atm_dct, bnd_dct)
    return xgr


def subgraph(xgr, atm_keys):
    """ return the subgraph induced by a subset of the atoms
    """
    assert set(atm_keys) <= set(atom_keys(xgr))
    bnd_keys = [bnd_key for bnd_key in bond_keys(xgr)
                if set(bnd_key) <= set(atm_keys)]
    atm_dct = _by_key(atoms(xgr), atm_keys)
    bnd_dct = _by_key(bonds(xgr), bnd_keys)
    sub_xgr = (atm_dct, bnd_dct)
    return sub_xgr


def delete_atoms(xgr, atm_keys):
    """ delete atoms from the molecular graph
    """
    all_atm_keys = set(atom_keys(xgr))
    atm_keys = set(atm_keys)
    assert atm_keys <= all_atm_keys
    atm_keys_left = all_atm_keys - atm_keys
    return subgraph(xgr, atm_keys_left)


def ring_keys_list(xgr):
    """ a series of key-sets for each ring in the graph
    """
    nxg = _nxg_from_graph(xgr)
    rng_keys_lst = _nxg_ring_keys_list(nxg)
    return rng_keys_lst


def backbone_keys(xgr):
    """ backbone keys of this molecular graph

    (backbone = heavy atoms and hydrogens that cannot be made implicit)
    """
    atm_keys = atom_keys(xgr)
    atm_sym_dct = atom_symbols(xgr)
    atm_ngb_keys_dct = atom_neighbor_keys(xgr)

    bbn_keys = []
    for atm_key in atm_keys:
        if atm_sym_dct[atm_key] != 'H':
            bbn_keys.append(atm_key)
        else:
            atm_ngb_keys = atm_ngb_keys_dct[atm_key]
            atm_ngb_syms = _values_by_key(atm_sym_dct, atm_ngb_keys)
            if all(sym == 'H' and atm_key < key
                   for key, sym in zip(atm_ngb_keys, atm_ngb_syms)):
                bbn_keys.append(atm_key)

    return tuple(bbn_keys)


def atom_explicit_hydrogen_keys(xgr):
    """ the number of explicit hydrogens, by atom
    """
    atm_sym_dct = atom_symbols(xgr)

    def _select_hydrogens(keys):
        return tuple(filter(lambda key: atm_sym_dct[key] == 'H', keys))

    atm_exp_hyd_keys = _transform_values(atom_neighbor_keys(xgr),
                                         func=_select_hydrogens)

    return atm_exp_hyd_keys


def hide_explicit_hydrogens(xgr, atm_keys=None):
    """ hide explicit hydrogens at these atoms (must be backbone atoms)
    """
    atm_keys = atm_keys if atm_keys is not None else backbone_keys(xgr)
    assert set(atm_keys) <= set(backbone_keys(xgr))

    atm_exp_hyd_keys_lst = _values_by_key(
        atom_explicit_hydrogen_keys(xgr), atm_keys)

    atm_imp_hyd_cnts = _values_by_key(atom_hydrogen_counts(xgr), atm_keys)
    atm_exp_hyd_cnts = tuple(map(len, atm_exp_hyd_keys_lst))
    atm_tot_hyd_cnt_dct = dict(
        zip(atm_keys, numpy.add(atm_imp_hyd_cnts, atm_exp_hyd_cnts)))
    xgr = change_atom_hydrogen_counts(xgr, atm_tot_hyd_cnt_dct)

    exp_hyd_keys = tuple(_chain(*atm_exp_hyd_keys_lst))
    xgr = delete_atoms(xgr, exp_hyd_keys)
    return xgr


def connectivity_graph(xgr):
    """ connectivity graph of this molecular graph
    """
    return _from_data(
        atm_keys=atom_keys(xgr),
        bnd_keys=bond_keys(xgr),
        atm_dcts=[atom_symbols(xgr), atom_hydrogen_counts(xgr)],
        bnd_dct={}
    )


def highspin_resonance_graph(xgr):
    """ high-spin (all sigma bond) resonance graph of this molecular graph
    """
    bnd_keys = bond_keys(xgr)
    return _from_data(
        atm_keys=atom_keys(xgr),
        bnd_keys=bnd_keys,
        atm_dcts=[atom_symbols(xgr), atom_hydrogen_counts(xgr)],
        bnd_dct=_by_key({}, bnd_keys, fill=1)
    )


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
    atm_dcts = [_by_key(atm_dct, atm_keys) for atm_dct in atm_dcts]
    bnd_dct = _by_key(bnd_dct, bnd_keys)

    assert all(set(atm_dct.keys()) == set(atm_keys) for atm_dct in atm_dcts)
    assert set(bnd_dct.keys()) == set(bnd_keys)

    atm_dct = dict(zip(
        atm_keys,
        zip(*(_values_by_key(atm_dct, atm_keys) for atm_dct in atm_dcts))
    ))

    xgr = (atm_dct, bnd_dct)
    return xgr