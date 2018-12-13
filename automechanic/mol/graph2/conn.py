""" connectivity graph library; by analogy to InChI connectivity layer

vertices: atomic symbols and hydrogen counts
    (('O', 1), ('C', 2), ('C', 2), ...)
edges: bond connectivity only (no bond orders or anything)
    {{0, 1}: None, {1, 2}: None, ...}
"""
import numpy
from .base import vertices as _vertices
from .base import edges as _edges
from .base import vertex_keys as _vertex_keys
from .base import vertex_edges as _vertex_edges
from ..atom import valence as _atom_valence
from ..molfile import FMT as _MLF_FMT


def valences(cgr):
    """ valences for each atom (independent of connectivity)
    """
    atm_syms, _ = zip(*_vertices(cgr))
    atm_val_elec_cnts = tuple(map(_atom_valence, atm_syms))
    return atm_val_elec_cnts


def sigma_bond_counts(cgr):
    """ sigma bond count for each atom
    """
    atm_keys = _vertex_keys(cgr)
    _, atm_hcnts = zip(*_vertices(cgr))
    assert all(numpy.greater_equal(atm_hcnts, 0))
    atm_cnn_elec_cnts = tuple(len(_vertex_edges(cgr, key)) for key in atm_keys)
    atm_sig_elec_cnts = numpy.add(atm_cnn_elec_cnts, atm_hcnts)
    return tuple(atm_sig_elec_cnts)


def pi_or_radical_electron_counts(cgr):
    """ the number of pi or radical electrons for each atom
    """
    atm_pi_or_rad_elec_cnts = numpy.subtract(valences(cgr),
                                             sigma_bond_counts(cgr))
    # disallow hypervalent atoms
    assert all(numpy.greater_equal(atm_pi_or_rad_elec_cnts, 0))
    return tuple(atm_pi_or_rad_elec_cnts)


def possible_spin_multiplicities(cgr):
    """ possible spin multiplicities for this molecule connectivity
    """
    pi_or_rad_elec_cnt = sum(pi_or_radical_electron_counts(cgr))
    return _spin_multiplicities(pi_or_rad_elec_cnt)


def molfile(cgr):
    """ MOLFile string of a connectivity graph
    """
    atms = _vertices(cgr)
    cnns = _edges(cgr)
    atm_keys = _vertex_keys(cgr)
    atm_sig_cnts = sigma_bond_counts(cgr)
    atm_rad_cnts = pi_or_radical_electron_counts(cgr)
    atm_cnt = len(atms)
    bnd_cnt = len(cnns)
    is_chiral = False
    counts = _MLF_FMT.COUNTS.LINE(
        **{_MLF_FMT.COUNTS.NA_KEY: atm_cnt,
           _MLF_FMT.COUNTS.NB_KEY: bnd_cnt,
           _MLF_FMT.COUNTS.CHI_KEY: is_chiral})
    atom = ''
    for atm_key, (atm_sym, _), atm_sig_cnt, atm_rad_cnt in zip(
            atm_keys, atms, atm_sig_cnts, atm_rad_cnts):
        atom += _MLF_FMT.ATOM.LINE(**{_MLF_FMT.ATOM.I_KEY: atm_key + 1,
                                      _MLF_FMT.ATOM.S_KEY: atm_sym,
                                      _MLF_FMT.ATOM.X_KEY: 0,
                                      _MLF_FMT.ATOM.Y_KEY: 0,
                                      _MLF_FMT.ATOM.Z_KEY: 0,
                                      _MLF_FMT.ATOM.VAL_KEY: atm_sig_cnt,
                                      _MLF_FMT.ATOM.MULT_KEY: atm_rad_cnt+1})
    # TODO: 2. hydrogens
    # TODO: 3. bonds
    mlf = _MLF_FMT.STRING(**{_MLF_FMT.COUNTS_KEY: counts,
                             _MLF_FMT.ATOM_KEY: atom,
                             _MLF_FMT.BOND_KEY: 'bond\n'})
    print(mlf)
    raise NotImplementedError


def _spin_multiplicities(nelecs):
    assert nelecs >= 0
    mult_max = nelecs + 1
    mult_min = 2 if (mult_max % 2 == 0) else 1
    mults = tuple(range(mult_min, mult_max+1, 2))
    return mults
