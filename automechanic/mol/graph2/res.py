""" resonance graph module

atm_dct: {atm_key: (atm_sym, atm_hyd_cnt), ...}
bnd_dct: {bnd_key: bnd_cnt, ...}
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
from ._shared import highspin_resonance_graph
from ._shared import _from_data
from ._dict import values_by_key as _values_by_key
from ._molfile import from_data as _mlf_from_data
from ._irdkit import from_molfile as _rdm_from_molfile
from ._irdkit import to_inchi as _rdm_to_inchi
from .._inchi_aux import numbering as _ich_aux_numbering


def from_data(atm_keys, bnd_keys, atm_sym_dct, atm_hyd_cnt_dct, bnd_ord_dct):
    """ resonance graph from data
    """
    assert len(atm_keys) == len(atm_sym_dct), len(atm_hyd_cnt_dct)
    assert len(bnd_keys) == len(bnd_ord_dct)
    return _from_data(
        atm_keys=atm_keys,
        bnd_keys=bnd_keys,
        atm_dcts=[atm_sym_dct, atm_hyd_cnt_dct],
        bnd_dct=bnd_ord_dct
    )


def bond_orders(rgr):
    """ bond orders (alias of bonds)
    """
    return bonds(rgr)


def atom_bond_valences(rgr):
    """ bond valences, by atom
    """
    atm_keys = atom_keys(rgr)
    atm_hyd_cnts = _values_by_key(atom_hydrogen_counts(rgr), atm_keys)
    atm_bnd_dcts = _values_by_key(atom_bonds(rgr), atm_keys)
    atm_bnd_cnts = [sum(atm_bnd_dct.values()) for atm_bnd_dct in atm_bnd_dcts]
    atm_bnd_vlcs = numpy.add(atm_hyd_cnts, atm_bnd_cnts)
    return dict(zip(atm_keys, atm_bnd_vlcs))


def atom_radical_valences(rgr):
    """ radical valences, by atom
    """
    atm_keys = atom_keys(rgr)
    atm_bnd_vlcs = _values_by_key(atom_bond_valences(rgr), atm_keys)
    atm_tot_vlcs = _values_by_key(atom_total_valences(rgr), atm_keys)
    atm_rad_vlcs = numpy.subtract(atm_tot_vlcs, atm_bnd_vlcs)
    return {atm_key: atm_rad_vlc
            for atm_key, atm_rad_vlc in zip(atm_keys, atm_rad_vlcs)
            if atm_rad_vlc}


def inchi(rgr):
    """ InChI string of a resonance graph """
    ich, _ = _inchi_with_atom_priorities(rgr)
    return ich


def atom_inchi_numbers(rgr):
    """ InChI numbering of backbone atoms
    """
    _, nums = _inchi_with_atom_priorities(rgr)
    return nums


def _inchi_with_atom_priorities(rgr):
    """ InChI string of this resonance graph
    """
    atm_keys = atom_keys(rgr)
    bnd_keys = bond_keys(rgr)
    atm_syms = _values_by_key(atom_symbols(rgr), atm_keys)
    atm_bnd_vlcs = _values_by_key(atom_bond_valences(rgr), atm_keys)
    atm_rad_vlcs = _values_by_key(atom_radical_valences(rgr), atm_keys, fill=0)
    bnd_ords = _values_by_key(bond_orders(rgr), bnd_keys)
    mlf = _mlf_from_data(atm_keys, bnd_keys, atm_syms, atm_bnd_vlcs,
                         atm_rad_vlcs, bnd_ords)
    rdm = _rdm_from_molfile(mlf, with_stereo=False)
    ich, ich_aux = _rdm_to_inchi(rdm, with_aux_info=True)
    nums = _ich_aux_numbering(ich_aux)
    assert len(atm_keys) == len(nums)
    atm_ich_num_dct = dict(zip(atm_keys, nums))
    return ich, atm_ich_num_dct


__all__ = ['atoms', 'bonds', 'atom_keys', 'bond_keys', 'atom_symbols',
           'atom_hydrogen_counts', 'atom_nuclear_charges',
           'atom_total_valences', 'atom_bonds', 'atom_neighbor_keys',
           'branch', 'relabel', 'subgraph',
           'isomorphism', 'isomorphic', 'highspin_resonance_graph']
