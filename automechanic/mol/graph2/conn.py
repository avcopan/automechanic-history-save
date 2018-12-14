""" connectivity graph library; by analogy to InChI connectivity layer

vertices: atomic symbols and implicit hydrogen counts
    (('O', 1), ('C', 2), ('C', 2), ...)
edges: bond connectivity only (no bond orders or anything)
    {{0, 1}: None, {1, 2}: None, ...}
"""
import numpy
from .base import vertices as _vertices
from .base import vertex_keys as _vertex_keys
from .base import edge_keys as _edge_keys
from .base import vertex_edges as _vertex_edges
from .base import vertex_neighbor_keys as _vertex_neighbor_keys
from .base import permute_vertices as _permute_vertices
from ..atom import valence as _atom_valence
from ..molfile import FMT as _MLF


def atomic_symbols(cgr):
    """ atomic symbols
    """
    atm_syms, _ = zip(*_vertices(cgr))
    return atm_syms


def implicit_hydrogen_counts(cgr):
    """ atomic symbols
    """
    _, atm_hcnts = zip(*_vertices(cgr))
    return atm_hcnts


def valences(cgr):
    """ valences for each atom (independent of connectivity)
    """
    atm_syms = atomic_symbols(cgr)
    atm_val_elec_cnts = tuple(map(_atom_valence, atm_syms))
    return atm_val_elec_cnts


def _backbone_atom_keys(cgr):
    atm_keys = _vertex_keys(cgr)
    atm_syms = atomic_symbols(cgr)
    hv_atm_keys = tuple(key for key in atm_keys if atm_syms[key] != 'H')
    hy_atm_keys = tuple(key for key in atm_keys if atm_syms[key] == 'H'
                        if all(key < nkey and atm_syms[nkey] == 'H'
                               for nkey in _vertex_neighbor_keys(cgr, key)))
    return hv_atm_keys + hy_atm_keys


def make_implicit(cgr):
    """ make explicit hydrogens implicit
    """
    atm_keys_permutation = tuple(numpy.argsort(
        # should use numpy.equal, but it's not implemented for strings -- ??
        numpy.array(atomic_symbols(cgr)) == 'H'))
    cgr = _permute_vertices(cgr, atm_keys_permutation)
    atm_keys = _backbone_atom_keys(cgr)

    # atm_keys = _vertex_keys(cgr)
    # atm_syms = atomic_symbols(cgr)
    # atm_hcnts = implicit_hydrogen_counts(cgr)
    # print(atm_syms)
    # print(atm_hcnts)


def sigma_bond_counts(cgr):
    """ sigma bond count for each atom
    """
    atm_keys = _vertex_keys(cgr)
    atm_hcnts = implicit_hydrogen_counts(cgr)
    assert all(numpy.greater_equal(atm_hcnts, 0))
    atm_cnn_elec_cnts = tuple(len(_vertex_edges(cgr, key)) for key in atm_keys)
    atm_sig_elec_cnts = numpy.add(atm_cnn_elec_cnts, atm_hcnts)
    return tuple(atm_sig_elec_cnts)


def nonsigma_electron_counts(cgr):
    """ the number of pi or radical electrons for each atom
    """
    atm_nonsig_elec_cnts = numpy.subtract(valences(cgr),
                                          sigma_bond_counts(cgr))
    # disallow hypervalent atoms
    assert all(numpy.greater_equal(atm_nonsig_elec_cnts, 0))
    return tuple(atm_nonsig_elec_cnts)


def pi_or_radical_atoms(cgr):
    """ atoms with pi or radical electrons
    """
    return {atm_key: cnt for atm_key, cnt in
            zip(_vertex_keys(cgr), nonsigma_electron_counts(cgr)) if cnt > 0}


def possible_spin_multiplicities(cgr):
    """ possible spin multiplicities for this molecule connectivity
    """
    nonsig_elec_cnt = sum(nonsigma_electron_counts(cgr))
    return _spin_multiplicities(nonsig_elec_cnt)


def molfile(cgr):
    """ MOLFile string of a connectivity graph
    """
    atm_keys = _vertex_keys(cgr)
    bnd_keys = _edge_keys(cgr)
    atm_syms = atomic_symbols(cgr)
    atm_rad_cnts = nonsigma_electron_counts(cgr)
    atm_sig_cnts = sigma_bond_counts(cgr)

    # counts line
    atm_cnt = len(atm_keys)
    bnd_cnt = len(bnd_keys)
    is_chiral = False
    counts_line = _MLF.COUNTS.LINE(
        **{_MLF.COUNTS.NA_KEY: atm_cnt,
           _MLF.COUNTS.NB_KEY: bnd_cnt,
           _MLF.COUNTS.CHI_KEY: is_chiral})

    atom_block = ''.join((
        _MLF.ATOM.LINE(**{_MLF.ATOM.I_KEY: atm_key+1,
                          _MLF.ATOM.S_KEY: atm_sym,
                          _MLF.ATOM.X_KEY: 0,
                          _MLF.ATOM.Y_KEY: 0,
                          _MLF.ATOM.Z_KEY: 0,
                          _MLF.ATOM.VAL_KEY: atm_sig_cnt,
                          _MLF.ATOM.MULT_KEY: atm_rad_cnt+1,
                          _MLF.ATOM.CFG_KEY: 0})
        for atm_key, atm_sym, atm_rad_cnt, atm_sig_cnt
        in zip(atm_keys, atm_syms, atm_rad_cnts, atm_sig_cnts)))

    bond_block = ''.join((
        _MLF.BOND.LINE(**{_MLF.BOND.I_KEY: i+1,
                          _MLF.BOND.ORDER_KEY: 1,
                          _MLF.BOND.I1_KEY: min(bnd_key)+1,
                          _MLF.BOND.I2_KEY: max(bnd_key)+1})
        for i, bnd_key in enumerate(bnd_keys)))

    mlf = _MLF.STRING(**{_MLF.COUNTS_KEY: counts_line,
                         _MLF.ATOM_KEY: atom_block,
                         _MLF.BOND_KEY: bond_block})
    return mlf


def _spin_multiplicities(nelecs):
    assert nelecs >= 0
    mult_max = nelecs + 1
    mult_min = 2 if (mult_max % 2 == 0) else 1
    mults = tuple(range(mult_min, mult_max+1, 2))
    return mults
