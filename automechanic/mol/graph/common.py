""" functions common to all graph types

terminology:
    - vertex key: an integer list index
    - edge key: frozenset of two vertex keys
    - vertices: tuple, keyed by vertex key
    - edges: dict, keyed by edge key
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
