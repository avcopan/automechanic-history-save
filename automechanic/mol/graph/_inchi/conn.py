""" connectivity graph -> InChI string conversion
"""
from .molfile import from_data as _mlf_from_data
from .inchi_aux import sorted_atom_keys as _ich_aux_sorted_atom_keys
from ._rdkit import from_molfile as _rdm_from_molfile
from ._rdkit import to_inchi_with_aux_info as _rdm_to_inchi_with_aux_info
from .._base import atom_keys
from .._base import bond_keys
from .._base import atom_symbols
from .._base import atom_bond_valences
from .._base import atom_radical_valences
from .._base import backbone_keys
from .._base import atom_explicit_hydrogen_keys
from .._base import bond_orders
from .._dict import values_by_key as _values_by_key


def inchi(xgr):
    """ InChI string of this connectivity graph
    """
    ich, _ = _inchi_with_numbering(xgr)
    return ich


def atom_inchi_numbers(xgr):
    """ InChI numbers, by atom
    """
    _, atm_ich_num_dct = _inchi_with_numbering(xgr)
    return atm_ich_num_dct


def _inchi_with_numbering(xgr):
    atm_keys = atom_keys(xgr)
    bnd_keys = bond_keys(xgr)
    atm_syms = _values_by_key(atom_symbols(xgr), atm_keys)
    atm_bnd_vlcs = _values_by_key(atom_bond_valences(xgr), atm_keys)
    atm_rad_vlcs = _values_by_key(atom_radical_valences(xgr), atm_keys)
    bnd_ords = _values_by_key(bond_orders(xgr), bnd_keys)
    mlf, mlf_atm_key_dct = _mlf_from_data(atm_keys, bnd_keys, atm_syms,
                                          atm_bnd_vlcs, atm_rad_vlcs, bnd_ords)
    rdm = _rdm_from_molfile(mlf)
    ich, ich_aux = _rdm_to_inchi_with_aux_info(rdm)

    # determine the inchi numbering from the AuxInfo string
    ich_srt_mlf_bbn_keys = _ich_aux_sorted_atom_keys(ich_aux)
    ich_srt_bbn_keys = _values_by_key(mlf_atm_key_dct, ich_srt_mlf_bbn_keys)
    assert set(ich_srt_bbn_keys) == set(backbone_keys(xgr))
    atm_exp_hyd_keys_dct = atom_explicit_hydrogen_keys(xgr)
    last_num = len(ich_srt_bbn_keys)
    atm_ich_num_dct = {}
    for num, bbn_key in enumerate(ich_srt_bbn_keys):
        atm_ich_num_dct[bbn_key] = num
        for exp_hyd_key in atm_exp_hyd_keys_dct[bbn_key]:
            atm_ich_num_dct[exp_hyd_key] = last_num
            last_num += 1
    assert set(atm_ich_num_dct.keys()) == set(atm_keys)

    return ich, atm_ich_num_dct
