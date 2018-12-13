""" base graph library; functions common to all graph types

terminology:
    - vertices: a tuple; tuple indices are called `vertex keys`
    - edges: a dictionary; keyed by frozenset pairs of `vertex keys`, which
      will be termed `edge keys`
"""


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
