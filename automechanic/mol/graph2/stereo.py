""" stereo graph module

atm_dct: {atm_key: (atm_sym, atm_hyd_cnt, atm_par), ...}
bnd_dct: {bnd_key: bnd_par, ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
from ._shared import atoms
from ._shared import bonds
from ._shared import atom_keys
from ._shared import bond_keys
from ._shared import atom_symbols
from ._shared import atom_hydrogen_counts
from ._shared import atom_neighbor_keys
from ._shared import relabel
from ._shared import isomorphism
from ._shared import isomorphic
from ._shared import _atom_stereo_parities
from ._shared import _from_data
from ._molfile import FMT as _MLF_FMT


def from_data(atm_keys, bnd_keys, atm_sym_dct, atm_hyd_cnt_dct, atm_par_dct,
              bnd_par_dct):
    """ resonance graph from data
    """
    assert len(atm_keys) == len(atm_sym_dct), len(atm_hyd_cnt_dct)
    return _from_data(
        atm_keys=atm_keys,
        bnd_keys=bnd_keys,
        atm_dcts=[atm_sym_dct, atm_hyd_cnt_dct, atm_par_dct],
        bnd_dct=bnd_par_dct
    )


def atom_parities(sgr):
    """ atom stereo parities
    """
    return {atm_key: atm_par
            for atm_key, atm_par in _atom_stereo_parities(sgr).items()
            if atm_par is not None}


def bond_parities(sgr):
    """ bond stereo parities
    """
    return {bnd_key: bnd_par
            for bnd_key, bnd_par in bonds(sgr).items()
            if bnd_par is not None}


def reflection(sgr):
    """ stereo graph reflection (inverts atom parities)
    """
    return from_data(
        atm_keys=atom_keys(sgr),
        bnd_keys=bond_keys(sgr),
        atm_sym_dct=atom_symbols(sgr),
        atm_hyd_cnt_dct=atom_hydrogen_counts(sgr),
        atm_par_dct={atm_key: not atm_par
                     for atm_key, atm_par in atom_parities(sgr).items()},
        bnd_par_dct=bond_parities(sgr)
    )


def is_chiral(sgr):
    """ is this stereo graph chiral?
    """
    return isomorphic(sgr, reflection(sgr))


def inchi(sgr):
    """ InChI string of this stereo graph
    """
    atm_keys = atom_keys(sgr)
    bnd_keys = bond_keys(sgr)
    is_chi = is_chiral(sgr)
    bnd_cfgs = _molfile_bond_configurations(sgr)
    print(atm_keys)
    print(bnd_keys)
    print(bnd_cfgs)
    print(is_chi)
    print(_MLF_FMT)


def _molfile_bond_configurations(sgr):
    atm_keys, atm_pars = zip(*atom_parities(sgr).items())
    atm_ngb_keys_dct = atom_neighbor_keys(sgr)
    atm_ngb_keys = tuple(map(next, map(iter, map(sorted, map(
        atm_ngb_keys_dct.__getitem__, atm_keys)))))
    return {frozenset({atm_key, atm_ngb_key}): atm_par
            for atm_key, atm_ngb_key, atm_par
            in zip(atm_keys, atm_ngb_keys, atm_pars)}


__all__ = ['atoms', 'bonds', 'atom_keys', 'bond_keys', 'atom_symbols',
           'atom_hydrogen_counts', 'atom_neighbor_keys', 'relabel',
           'isomorphism', 'isomorphic']
