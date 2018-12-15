""" functions operating on connectivity graphs

terminology:
    - atm_syms: atomic symbols, e.g. ('C', 'H', 'H',...)
    - bnd_cnns: bond connections only, e.g. {{0, 1}: None, {0, 2}: None, ...}
"""
from .common import vertices
from .common import vertex_keys
from .common import vertex_edges
from ..atom import valence as _atom_valence


def possible_spin_multiplicities(cgr):
    """ possible spin multiplicities for this connectivity graph
    """
    atm_syms = vertices(cgr)
    atm_keys = vertex_keys(cgr)
    nsigma_elecs = sum(len(vertex_edges(cgr, atm_key)) for atm_key in atm_keys)
    nvalence_elecs = sum(_atom_valence(atm_sym) for atm_sym in atm_syms)
    return _spin_multiplicities(nvalence_elecs - nsigma_elecs)


def _spin_multiplicities(nelecs):
    mult_max = nelecs + 1
    mult_min = 2 if (mult_max % 2 == 0) else 1
    mults = tuple(range(mult_min, mult_max+1, 2))
    return mults
