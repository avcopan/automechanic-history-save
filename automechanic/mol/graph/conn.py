""" connectivity graph library; by analogy to InChI connectivity layer

vertices: atomic symbols and implicit hydrogen counts
    (('O', 1), ('C', 2), ('C', 2), ...)
edges: bond connectivity only (no bond orders or anything)
    {{0, 1}: None, {1, 2}: None, ...}
"""
from itertools import starmap as _starmap
from itertools import combinations as _combinations
from itertools import permutations as _permutations
from .base import vertices as _vertices
from .base import edge_keys as _edge_keys
from .base import vertex_neighbor_keys as _vertex_neighbor_keys
from .base import branch as _branch
from .base import cycle_keys_list as _cycle_keys_list
from .base import isomorphic as _isomorphic
from ._shared import atomic_symbols
from ._shared import hydrogen_counts
from ._shared import valences
from ._shared import change_hydrogen_count
from ._shared import backbone_keys
from ._shared import nonbackbone_keys
from ._shared import neighboring_hydrogen_keys
from ._shared import make_hydrogens_implicit
from .res import radical_electron_counts as _res_radical_electron_counts
from .res import sigma_bond_counts as _res_sigma_bond_counts
from .res import (possible_spin_multiplicities as
                  _res_possible_spin_multiplicities)
from .res import inchi as _res_inchi
from .res import inchi_numbering as _res_inchi_numbering
from .res import open_pi_bond_keys as _res_open_pi_bond_keys
from .res import lowspin_resonances as _res_lowspin_resonances


def _no_pi_resonance_graph(cgr):
    atms = _vertices(cgr)
    bnds = {bnd_key: 1 for bnd_key in _edge_keys(cgr)}
    rgr = (atms, bnds)
    return rgr


def sigma_bond_counts(cgr):
    """ sigma bond count for each atom
    """
    rgr = _no_pi_resonance_graph(cgr)
    return _res_sigma_bond_counts(rgr)


def nonsigma_electron_counts(cgr):
    """ the number of pi or radical electrons for each atom
    """
    rgr = _no_pi_resonance_graph(cgr)
    return _res_radical_electron_counts(rgr)


def potential_pi_bond_keys(cgr):
    """ keys for connections that could be pi bonds
    """
    rgr = _no_pi_resonance_graph(cgr)
    return _res_open_pi_bond_keys(rgr)


def lowspin_resonance_graphs(cgr):
    """ all possible resonance graphs with this sigma bonding structure
    """
    rgr = _no_pi_resonance_graph(cgr)
    return _res_lowspin_resonances(rgr)


def resonance_graph(cgr):
    """ an arbitrary low-spin resonance graph with this sigma bonding structure
    """
    rgrs = lowspin_resonance_graphs(cgr)
    assert rgrs
    return rgrs[0]


def possible_spin_multiplicities(cgr):
    """ possible spin multiplicities for this molecule connectivity
    """
    rgr = _no_pi_resonance_graph(cgr)
    return _res_possible_spin_multiplicities(rgr)


def inchi(cgr):
    """ InChI string of a connectivity graph
    """
    return _res_inchi(resonance_graph(cgr))


def inchi_numbering(cgr):
    """ InChI numbering of backbone atoms
    """
    return _res_inchi_numbering(resonance_graph(cgr))


def stereogenic_atom_keys(cgr):
    """ stereogenic atoms and their neighbors in CIP priority
    """
    print(cgr)


def stereogenic_atoms(cgr):
    """ identify all stereogenic atoms in the connectivity graph

    currently missing allene stereo
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

    currently missing cumulene stereo
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

    ret = tuple(bnd_key for bnd_key in filter(_bond_is_stereo_candidate,
                                              potential_pi_bond_keys(cgr))
                if not any(_starmap(_atom_is_symmetric_on_bond,
                                    _permutations(bnd_key))))
    return ret


__all__ = [
    'atomic_symbols', 'hydrogen_counts', 'valences', 'change_hydrogen_count',
    'backbone_keys', 'nonbackbone_keys', 'neighboring_hydrogen_keys',
    'make_hydrogens_implicit'
]
