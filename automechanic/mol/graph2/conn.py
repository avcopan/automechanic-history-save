""" connectivity graph library; by analogy to InChI connectivity layer

vertices: atomic symbols and implicit hydrogen counts
    (('O', 1), ('C', 2), ('C', 2), ...)
edges: bond connectivity only (no bond orders or anything)
    {{0, 1}: None, {1, 2}: None, ...}
"""
from itertools import starmap as _starmap
from itertools import combinations as _combinations
from itertools import permutations as _permutations
import numpy
from .base import vertices as _vertices
from .base import edges as _edges
from .base import vertex_keys as _vertex_keys
from .base import edge_keys as _edge_keys
from .base import vertex_edges as _vertex_edges
from .base import vertex_neighbor_keys as _vertex_neighbor_keys
from .base import delete_vertices as _delete_vertices
from .base import branch as _branch
from .base import cycle_keys_list as _cycle_keys_list
from .base import permute_vertices as _permute_vertices
from .base import isomorphic as _isomorphic
from ..atom import valence as _atom_valence
from ..atom import lone_pair_count as _atom_lone_pair_count
from ..molfile import FMT as _MLF


def atomic_symbols(cgr):
    """ atomic symbols
    """
    syms, _ = zip(*_vertices(cgr))
    return syms


def hydrogen_counts(cgr):
    """ atomic symbols
    """
    _, hcnts = zip(*_vertices(cgr))
    return hcnts


def valences(cgr):
    """ valences for each atom (independent of connectivity)
    """
    syms = atomic_symbols(cgr)
    val_elec_cnts = tuple(map(_atom_valence, syms))
    return val_elec_cnts


def lone_pair_counts(cgr):
    """ lone pair counts for each atom (independent of connectivity)
    """
    syms = atomic_symbols(cgr)
    lp_cnts = tuple(map(_atom_lone_pair_count, syms))
    return lp_cnts


def change_hydrogen_count(cgr, key, nhyd):
    """ change the hydrogen count of an atom
    """
    syms = atomic_symbols(cgr)
    hcnts = list(hydrogen_counts(cgr))
    edgs = _edges(cgr)
    hcnts[key] += nhyd
    atms = tuple(zip(syms, hcnts))
    return (atms, edgs)


def _is_backbone_key(cgr, key):
    assert key in _vertex_keys(cgr)
    syms = atomic_symbols(cgr)
    sym = syms[key]
    nkeys = _vertex_neighbor_keys(cgr, key)
    nsyms = [syms[nkey] for nkey in nkeys]
    return sym != 'H' or all(
        nsym == 'H' and key < nkey
        for nsym, nkey in zip(nsyms, nkeys))


def _backbone_keys(cgr):
    return tuple(key for key in _vertex_keys(cgr)
                 if _is_backbone_key(cgr, key))


def _non_backbone_keys(cgr):
    return tuple(key for key in _vertex_keys(cgr)
                 if not _is_backbone_key(cgr, key))


def _neighboring_hydrogen_keys(cgr, key):
    assert key in _vertex_keys(cgr)
    syms = atomic_symbols(cgr)
    return tuple(nkey for nkey in _vertex_neighbor_keys(cgr, key)
                 if syms[nkey] == 'H')


def make_hydrogens_implicit(cgr):
    """ make explicit hydrogens implicit
    """
    keys_perm = _backbone_keys(cgr) + _non_backbone_keys(cgr)
    cgr = _permute_vertices(cgr, keys_perm)
    bbn_keys = _backbone_keys(cgr)
    for bbn_key in bbn_keys:
        hyd_keys = _neighboring_hydrogen_keys(cgr, bbn_key)
        nhyd = len(hyd_keys)
        cgr = _delete_vertices(cgr, hyd_keys)
        cgr = change_hydrogen_count(cgr, bbn_key, nhyd)
    return cgr


def sigma_bond_counts(cgr):
    """ sigma bond count for each atom
    """
    keys = _vertex_keys(cgr)
    hcnts = hydrogen_counts(cgr)
    assert all(numpy.greater_equal(hcnts, 0))
    cnn_elec_cnts = tuple(len(_vertex_edges(cgr, key)) for key in keys)
    sig_elec_cnts = numpy.add(cnn_elec_cnts, hcnts)
    return tuple(sig_elec_cnts)


def nonsigma_electron_counts(cgr):
    """ the number of pi or radical electrons for each atom
    """
    nonsig_elec_cnts = numpy.subtract(valences(cgr), sigma_bond_counts(cgr))
    # disallow hypervalent atoms
    assert all(numpy.greater_equal(nonsig_elec_cnts, 0))
    return tuple(nonsig_elec_cnts)


def pi_or_radical_atom_keys(cgr):
    """ atoms with pi or radical electrons
    """
    return tuple(key for key, cnt in enumerate(nonsigma_electron_counts(cgr))
                 if cnt > 0)


def potential_pi_bond_keys(cgr):
    """ potential pi bond keys (adjacent sites with pi/radical electrons)
    """
    bnd_keys = _edge_keys(cgr)
    pi_rad_keys = set(pi_or_radical_atom_keys(cgr))
    return tuple(bnd_key for bnd_key in bnd_keys if bnd_key <= pi_rad_keys)


def stereogenic_atoms(cgr):
    """ identify all stereogenic atoms in the connectivity graph
    """
    keys = []
    sig_bnd_cnts = sigma_bond_counts(cgr)
    for key, sig_bnd_cnt in enumerate(sig_bnd_cnts):
        if sig_bnd_cnt == 4:
            nkeys = set(_vertex_neighbor_keys(cgr, key))
            if len(nkeys) >= 3:
                bchs = [_branch(cgr, key, nkey, excl_keys=nkeys-{nkey})
                        for nkey in nkeys]
                if not any(_starmap(_isomorphic, _combinations(bchs, r=2))):
                    keys.append(key)
    return tuple(keys)


def stereogenic_bonds(cgr):
    """ identify all stereogenic bonds in this connectivity graph
    """
    def _bond_is_stereo_candidate(bnd_key):
        nsig_elec_cnts = nonsigma_electron_counts(cgr)
        not_linear = any(nsig_elec_cnts[atm_key] < 2 for atm_key in bnd_key)
        cyc_keys_lst = _cycle_keys_list(cgr)
        # follow InChI rule, ignoring stereo bonds in rings of size < 8
        not_in_small_cycle = not any(bnd_key < set(cyc_keys) for cyc_keys in
                                     cyc_keys_lst if len(cyc_keys) < 8)
        return not_linear and not_in_small_cycle

    def _atom_is_symmetric_on_bond(atm_key, bonded_atm_key):
        nkeys = _vertex_neighbor_keys(cgr, atm_key)
        assert bonded_atm_key in nkeys
        nkeys = set(nkeys) - {bonded_atm_key}
        assert 0 <= len(nkeys) <= 2
        ret = not bool(nkeys)
        if len(nkeys) == 2:
            nkey1, nkey2 = nkeys
            bch1 = _branch(cgr, atm_key, nkey1, excl_keys=(bonded_atm_key,))
            bch2 = _branch(cgr, atm_key, nkey2, excl_keys=(bonded_atm_key,))
            ret = _isomorphic(bch1, bch2)
        return ret

    return tuple(bnd_key for bnd_key in filter(_bond_is_stereo_candidate,
                                               potential_pi_bond_keys(cgr))
                 if not any(_starmap(_atom_is_symmetric_on_bond,
                                     _permutations(bnd_key))))


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
