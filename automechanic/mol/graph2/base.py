""" base graph library; functions common to all graph types

terminology:
    - vertices: a tuple; tuple indices are called `vertex keys`
    - edges: a dictionary; keyed by frozenset pairs of `vertex keys`, which
      will be termed `edge keys`
"""
from itertools import chain as _chain
from ._perm import inverse as _inverse_permutation


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


def vertex_edges(gra, vkey):
    """ edges at a given vertex
    """
    edgs = edges(gra)
    return {ekey: evals for ekey, evals in edgs.items() if vkey in ekey}


def vertex_neighbor_keys(gra, vkey):
    """ keys for vertices adjacent to this one
    """
    vtx_ekeys = vertex_edges(gra, vkey).keys()
    vtx_nkeys, = zip(*(ekey - {vkey} for ekey in vtx_ekeys))
    return tuple(sorted(vtx_nkeys))


def induced_subgraph(gra, vkeys):
    """ the subgraph induced by a set of vertices
    """
    assert set(vkeys) <= set(vertex_keys(gra))
    vtcs = vertices(gra)
    edgs = edges(gra)
    ekeys = tuple(ekey for ekey in edge_keys(gra)
                  if all(edg_vkey in vkeys for edg_vkey in ekey))
    vkey_map = dict(map(reversed, enumerate(vkeys)))
    vtcs = tuple(vtcs[vkey] for vkey in vkeys)
    edgs = {ekey: edgs[ekey] for ekey in ekeys}
    edgs = _edges_for_new_vertex_keys(edgs, vkey_map)
    return (vtcs, edgs)


def delete_vertices(gra, vkeys):
    """ delete vertices and their associated edges
    """
    all_vkeys = vertex_keys(gra)
    assert set(vkeys) <= set(all_vkeys)
    return induced_subgraph(gra, set(all_vkeys) - set(vkeys))


def _branch_vertex_keys(gra, vkeys, seen_vkeys, excl_vkeys):
    new_vkeys = vkeys - seen_vkeys
    if new_vkeys:
        new_vkey_neighbors = set(_chain(
            *(vertex_neighbor_keys(gra, vkey) for vkey in new_vkeys)))
        vkeys.update(new_vkey_neighbors - excl_vkeys)
        seen_vkeys.update(new_vkeys)
        _branch_vertex_keys(gra, vkeys, seen_vkeys, excl_vkeys)
    return vkeys


def branch(gra, vkey1, vkey2, excl_vkeys=None):
    """ return the branch extending from vkey1 and vkey2 away from vkey1

    vkey1 and vkey2 must be adjacent
    """
    assert frozenset({vkey1, vkey2}) in edge_keys(gra)
    excl_vkeys = {vkey1} if excl_vkeys is None else excl_vkeys.union({vkey1})
    vkeys = _branch_vertex_keys(gra, {vkey2}, set(), excl_vkeys)
    vkeys.add(vkey1)
    return induced_subgraph(gra, vkeys)


def permute_vertices(gra, vkeys_permutation):
    """ permute graph vertices
    """
    assert tuple(sorted(vkeys_permutation)) == vertex_keys(gra)
    vtcs = vertices(gra)
    edgs = edges(gra)
    ret_vtcs = tuple(vtcs[i] for i in vkeys_permutation)

    vkey_map = dict(enumerate(_inverse_permutation(vkeys_permutation)))
    ret_edgs = _edges_for_new_vertex_keys(edgs, vkey_map=vkey_map)

    return (ret_vtcs, ret_edgs)


def _edges_for_new_vertex_keys(edgs, vkey_map):
    if not edgs:
        ret_edgs = {}
    else:
        ekeys = edgs.keys()
        evals = edgs.values()
        ret_ekeys = tuple(frozenset(map(vkey_map.__getitem__, ekey))
                          for ekey in ekeys)
        ret_edgs = dict(zip(ret_ekeys, evals))
    return ret_edgs
