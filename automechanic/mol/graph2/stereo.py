""" stereo graph module

atm_dct: {atm_key: (atm_sym, atm_hyd_cnt, atm_par), ...}
bnd_dct: {bnd_key: bnd_par, ...}
bnd_key := frozenset({atm1_key, atm2_key})
"""
import numpy
from ._shared import atoms
from ._shared import bonds
from ._shared import atom_keys
from ._shared import bond_keys
from ._shared import atom_symbols
from ._shared import atom_hydrogen_counts
from ._shared import atom_nuclear_charges
from ._shared import atom_total_valences
from ._shared import atom_bonds
from ._shared import atom_neighbor_keys
from ._shared import branch
from ._shared import relabel
from ._shared import subgraph
from ._shared import isomorphism
from ._shared import isomorphic
from ._shared import ring_keys_list
from ._shared import highspin_resonance_graph
from ._shared import _from_data
from ._shared import _atom_stereo_parities
from ._dict import filter_by_value as _filter_by_value
from ._to_inchi import atom_inchi_numbers


def from_data(atm_keys, bnd_keys, atm_sym_dct, atm_hyd_cnt_dct, atm_par_dct,
              bnd_par_dct):
    """ stereo graph from data
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
    return _filter_by_value(_atom_stereo_parities(sgr),
                            func=lambda val: val is not None)


def bond_parities(sgr):
    """ bond stereo parities
    """
    return _filter_by_value(bonds(sgr),
                            func=lambda val: val is not None)


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


def _coordinates(sgr):
    assert not ring_keys_list(sgr)  # currently assumes no rings
    atm_keys = atom_keys(sgr)
    atm_ngb_keys_dct = atom_neighbor_keys(sgr)

    # atm_ich_num_dct = atom_inchi_numbers(sgr)
    # atm_pars = atom_parities(sgr)

    disps = [(0., 1., 0.), (1., 0., 0.), (-1., 0., 0.), (0., -1., 0.)]

    atm_xyz_dct = {}

    atm_key = atm_keys[0]
    atm_xyz = numpy.array([0., 0., 0.])
    atm_xyz_dct[atm_key] = atm_xyz

    atm_ngb_keys = atm_ngb_keys_dct[atm_key]
    assert len(disps) > len(atm_ngb_keys)
    for atm_ngb_key, disp in zip(atm_ngb_keys, disps):
        atm_xyz_dct[atm_ngb_key] = numpy.add(atm_xyz, disp)

    return atm_xyz_dct


def inchi(sgr):
    """ InChI string of this stereo graph
    """
    atm_xyz_dct = _coordinates(sgr)
    print(atm_xyz_dct)


__all__ = ['atoms', 'bonds', 'atom_keys', 'bond_keys', 'atom_symbols',
           'atom_hydrogen_counts', 'atom_nuclear_charges',
           'atom_total_valences', 'atom_bonds', 'atom_neighbor_keys',
           'branch', 'relabel', 'subgraph',
           'isomorphism', 'isomorphic', 'highspin_resonance_graph',
           'atom_inchi_numbers']
