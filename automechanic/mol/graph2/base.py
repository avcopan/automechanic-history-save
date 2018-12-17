""" base graph library; functions common to all graph types

terminology:
    - vertices: a tuple; tuple indices are called `vertex keys`
    - edges: a dictionary; keyed by frozenset pairs of `vertex keys`, which
      will be termed `edge keys`
"""
from itertools import chain as _chain
from ._perm import inverse as _inverse_permutation
from ._inetworkx import from_automechanic_graph as _nxg_from_graph
from ._inetworkx import isomorphism as _nxg_isomorphism


def vertices(gra):
    """ vertices
    """
    vtcs, _ = gra
    return tuple(vtcs)


def edges(gra):
    """ edges
    """
    _, edgs = gra
    return dict(edgs)


def vertex_keys(gra):
    """ vertex keys
    """
    return tuple(range(len(vertices(gra))))


def edge_keys(gra):
    """ edge keys
    """
    return set(iter(edges(gra)))


def vertex_edges(gra, key):
    """ edges at a given vertex
    """
    edgs = edges(gra)
    return {ekey: evals for ekey, evals in edgs.items() if key in ekey}


def vertex_neighbor_keys(gra, key):
    """ keys for vertices adjacent to this one
    """
    vtx_ekeys = vertex_edges(gra, key).keys()
    vtx_nkeys, = zip(*(ekey - {key} for ekey in vtx_ekeys))
    return tuple(sorted(vtx_nkeys))


def induced_subgraph(gra, keys):
    """ the subgraph induced by a set of vertices
    """
    assert set(keys) <= set(vertex_keys(gra))
    vtcs = vertices(gra)
    edgs = edges(gra)
    ekeys = tuple(ekey for ekey in edge_keys(gra)
                  if all(edg_vkey in keys for edg_vkey in ekey))
    vkey_map = dict(map(reversed, enumerate(keys)))
    vtcs = tuple(vtcs[vkey] for vkey in keys)
    edgs = {ekey: edgs[ekey] for ekey in ekeys}
    edgs = _edges_for_new_vertex_keys(edgs, vkey_map)
    return (vtcs, edgs)


def delete_vertices(gra, keys):
    """ delete vertices and their associated edges
    """
    all_keys = vertex_keys(gra)
    assert set(keys) <= set(all_keys)
    return induced_subgraph(gra, set(all_keys) - set(keys))


def _branch_keys_recursive(gra, keys, seen_keys, excl_keys):
    new_keys = keys - seen_keys
    if new_keys:
        new_key_neighbors = set(_chain(*(vertex_neighbor_keys(gra, key)
                                         for key in new_keys)))
        keys.update(new_key_neighbors - excl_keys)
        seen_keys.update(new_keys)
        _branch_keys_recursive(gra, keys, seen_keys, excl_keys)
    return keys


def branch_keys(gra, key1, key2, excl_keys=None):
    """ return keys for the branch extending from two adjacent vertices
    """
    assert frozenset({key1, key2}) in edge_keys(gra)
    excl_keys = {key1} if excl_keys is None else set(excl_keys).union({key1})
    keys = _branch_keys_recursive(gra, {key2}, set(), excl_keys)
    keys.add(key1)
    return tuple(sorted(keys))


def branch(gra, key1, key2, excl_keys=None):
    """ return the branch extending from two adjacent vertices
    """
    bch_keys = branch_keys(gra, key1, key2, excl_keys=excl_keys)
    return induced_subgraph(gra, bch_keys)


def _ring_keys_list_recursive(gra, ring_keys_lst, chain_keys_lst):
    next_chain_keys_lst = []
    for chain_keys in chain_keys_lst:
        nkeys = set(vertex_neighbor_keys(gra, chain_keys[-1]))
        if chain_keys[0] in nkeys:
            ring_keys_lst.append(chain_keys)
        for nkey in nkeys - {chain_keys[0], chain_keys[-2]}:
            next_chain_keys_lst.append(chain_keys + [nkey])
        _ring_keys_list_recursive(gra, ring_keys_lst, next_chain_keys_lst)
    return ring_keys_lst


def ring_keys_list(gra):
    """ return keys for all rings (cycles) in the graph
    """
    rng_keys_lst = []
    for ekey in edge_keys(gra):
        if not any(set(rng_keys) > ekey for rng_keys in rng_keys_lst):
            rng_keys_lst += _ring_keys_list_recursive(gra, [], [list(ekey)])

    print(rng_keys_lst)


def permute_vertices(gra, perm_keys):
    """ permute graph vertices
    """
    assert tuple(sorted(perm_keys)) == vertex_keys(gra)
    vtcs = vertices(gra)
    edgs = edges(gra)
    ret_vtcs = tuple(vtcs[i] for i in perm_keys)

    key_map = dict(enumerate(_inverse_permutation(perm_keys)))
    ret_edgs = _edges_for_new_vertex_keys(edgs, key_map=key_map)

    return (ret_vtcs, ret_edgs)


def _edges_for_new_vertex_keys(edgs, key_map):
    if not edgs:
        ret_edgs = {}
    else:
        ekeys = edgs.keys()
        evals = edgs.values()
        ret_ekeys = tuple(frozenset(map(key_map.__getitem__, ekey))
                          for ekey in ekeys)
        ret_edgs = dict(zip(ret_ekeys, evals))
    return ret_edgs


def isomorphic(gra1, gra2):
    """ whether or not these graphs are isomorphic
    """
    return isomorphism(gra1, gra2) is not None


def isomorphism(gra1, gra2):
    """ graph isomorphism (permutation of `gra1` to produce `gra2`)
    """
    nxg1 = _nxg_from_graph(gra1)
    nxg2 = _nxg_from_graph(gra2)
    iso = _nxg_isomorphism(nxg1, nxg2)
    return iso
