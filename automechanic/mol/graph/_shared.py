""" shared functions for connectivity, resonance, and stereo graphs
"""
from .base import vertices as _vertices
from .base import vertex_keys as _vertex_keys
from .base import vertex_neighbor_keys as _vertex_neighbor_keys
from .base import edges as _edges
from .base import edge_keys as _edge_keys
from .base import delete_vertices as _delete_vertices
from .base import permute_vertices as _permute_vertices
from ..atom import valence as _atom_valence


def atomic_symbols(xgr):
    """ atomic symbols
    """
    syms = list(zip(*_vertices(xgr)))[0]
    return syms


def hydrogen_counts(xgr):
    """ atomic symbols
    """
    hcnts = list(zip(*_vertices(xgr)))[1]
    return hcnts


def connectivity_graph(xgr):
    """ strip this graph down to a basic connectivity graph
    """
    syms = atomic_symbols(xgr)
    hcnts = hydrogen_counts(xgr)
    atms = tuple(zip(syms, hcnts))
    cnn_keys = _edge_keys(xgr)
    cnns = {cnn_key: None for cnn_key in cnn_keys}
    cgr = (atms, cnns)
    return cgr


def valences(cgr):
    """ valences for each atom (independent of connectivity)
    """
    syms = atomic_symbols(cgr)
    val_elec_cnts = tuple(map(_atom_valence, syms))
    return val_elec_cnts


def change_hydrogen_count(cgr, key, nhyd):
    """ change the hydrogen count of an atom
    """
    vtcs = _vertices(cgr)
    edgs = _edges(cgr)
    vtcs_lst = list(map(list, vtcs))
    vtcs_lst[key][1] += nhyd
    vtcs = tuple(map(tuple, vtcs_lst))
    return (vtcs, edgs)


def backbone_keys(cgr):
    """ molecule backbone atom keys

    (heavy atoms and necessarily explicit hydrogens)
    """
    return tuple(key for key in _vertex_keys(cgr)
                 if _is_backbone_key(cgr, key))


def nonbackbone_keys(cgr):
    """ molecule non-backbone atom keys

    (hydrogen atoms that could be made implicit)
    """
    return tuple(key for key in _vertex_keys(cgr)
                 if not _is_backbone_key(cgr, key))


def neighboring_hydrogen_keys(cgr, key):
    """ keys for neighboring hydrogen atoms
    """
    assert key in _vertex_keys(cgr)
    syms = atomic_symbols(cgr)
    return tuple(nkey for nkey in _vertex_neighbor_keys(cgr, key)
                 if syms[nkey] == 'H')


def _is_backbone_key(cgr, key):
    assert key in _vertex_keys(cgr)
    syms = atomic_symbols(cgr)
    sym = syms[key]
    nkeys = _vertex_neighbor_keys(cgr, key)
    nsyms = [syms[nkey] for nkey in nkeys]
    return sym != 'H' or all(
        nsym == 'H' and key < nkey
        for nsym, nkey in zip(nsyms, nkeys))


def make_hydrogens_implicit(cgr):
    """ make explicit hydrogens implicit
    """
    keys_perm = backbone_keys(cgr) + nonbackbone_keys(cgr)
    cgr = _permute_vertices(cgr, keys_perm)
    bbn_keys = backbone_keys(cgr)
    for bbn_key in bbn_keys:
        hyd_keys = neighboring_hydrogen_keys(cgr, bbn_key)
        nhyd = len(hyd_keys)
        cgr = _delete_vertices(cgr, hyd_keys)
        cgr = change_hydrogen_count(cgr, bbn_key, nhyd)
    return cgr
