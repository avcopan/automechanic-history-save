""" stereo graph library; by analogy to the InChI stereo layer

vertices: atomic symbols, implicit hydrogen counts, stereo parity
    (('O', 1, None), ('C', 1, None), ('C', 1, False), ...)
edges: bond connectivity and stereo parity
    {{0, 1}: None, {1, 2}: True, ...}

stereo atom values:
    None = no stereo
    False = negative-parity stereo
    True = positive-parity stereo
stereo bond values:
    None = no stereo
    False = negative-parity stereo
    True = positive-parity stereo
"""
from .base import vertex_keys as _vertex_keys
from .base import vertices as _vertices
from .base import edges as _edges
from ._shared import connectivity_graph


def atom_stereo_values(sgr):
    """ atom stereo values (None if not stereogenic, True/False otherwise)
    """
    pars = list(zip(*_vertices(sgr)))[2]
    return pars


def invert_atom_stereo(sgr, key):
    """ invert the parity of a stereogenic atom
    """
    atms = _vertices(sgr)
    cnns = _edges(sgr)
    atms_lst = list(map(list, atms))
    val = atms_lst[key][2]
    atms_lst[key][2] = not val if val is not None else None
    atms = tuple(map(tuple, atms_lst))
    return (atms, cnns)


def reflect(sgr):
    """ reflect the graph through a mirror plane, inverting all atom stereo
    """
    raise NotImplementedError(sgr)


__all__ = ['connectivity_graph']
