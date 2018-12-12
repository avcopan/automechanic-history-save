""" functions operating on resonance graphs

terminology:
    - atm_syms: atomic symbols, e.g. ('C', 'H', 'H',...)
    - bnd_ords: bond connections only, e.g. {{0, 1}: None, {0, 2}: None, ...}
"""
import numpy
from .common import vertices
from .common import edges
from .common import vertex_keys
from .common import vertex_edges
from ..atom import valence as _atom_valence


def radical_sites(rgr):
    """ radical electrons by atom
    """
    atm_syms = vertices(rgr)
    atm_keys = vertex_keys(rgr)
    bnd_ords = edges(rgr)
    nval_elecs_by_atom = [_atom_valence(atm_sym) for atm_sym in atm_syms]
    nbnd_elecs_by_atom = [
        sum(bnd_ords[bnd_key] for bnd_key in vertex_edges(rgr, atm_key))
        for atm_key in atm_keys]
    nrad_elecs_by_atom = numpy.subtract(nval_elecs_by_atom, nbnd_elecs_by_atom)
    rad_dct = {atm_key: nrad_elecs for atm_key, nrad_elecs
               in zip(atm_keys, nrad_elecs_by_atom) if nrad_elecs > 0}
    return rad_dct
