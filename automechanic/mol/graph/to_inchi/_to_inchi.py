""" graph -> InChI string conversion
"""
from ._molfile import from_data as _mlf_from_data
from ._inchi_aux import sorted_atom_keys as _ich_aux_sorted_atom_keys
from ._rdkit import from_molfile as _rdm_from_molfile
from ._rdkit import to_inchi_with_aux_info as _rdm_to_inchi_with_aux_info
from .._base import atom_keys
from .._base import bond_keys
from .._base import atom_symbols
from .._base import implicit
from .._res import atom_bond_valences
from .._res import atom_radical_valences
from .._res import lowspin_resonance
from .._base import backbone_keys
from .._base import atom_explicit_hydrogen_keys
from .._base import bond_orders
from .._dict import values_by_key as _values_by_key


def with_numbering(cgr, is_chi=False, atm_xyz_dct=None):
    """ InChI string with numbering from a connectivity graph

    For stereo InChIs, pass cartesian coordinates and set the chirality flag.
    """
    ret = None
    if False:
        pass
    else:
        ret = _with_numbering(cgr, is_chi=False, atm_xyz_dct=None)

    return ret


def _exceptions(cgr):
    cgr = implicit(cgr)


def _with_numbering(cgr, is_chi=False, atm_xyz_dct=None):
    rgr = lowspin_resonance(cgr)
    atm_keys = atom_keys(rgr)
    bnd_keys = bond_keys(rgr)
    atm_syms = _values_by_key(atom_symbols(rgr), atm_keys)
    atm_bnd_vlcs = _values_by_key(atom_bond_valences(rgr), atm_keys)
    atm_rad_vlcs = _values_by_key(atom_radical_valences(rgr), atm_keys)
    atm_xyzs = (None if atm_xyz_dct is None else
                _values_by_key(atm_xyz_dct, atm_keys))
    bnd_ords = _values_by_key(bond_orders(rgr), bnd_keys)
    mlf, mlf_atm_key_dct = _mlf_from_data(atm_keys, bnd_keys, atm_syms,
                                          atm_bnd_vlcs, atm_rad_vlcs, bnd_ords,
                                          is_chi=is_chi, atm_xyzs=atm_xyzs)
    rdm = _rdm_from_molfile(mlf)
    ich, ich_aux = _rdm_to_inchi_with_aux_info(rdm)

    # determine the inchi numbering from the AuxInfo string
    ich_srt_mlf_bbn_keys = _ich_aux_sorted_atom_keys(ich_aux)
    ich_srt_bbn_keys = _values_by_key(mlf_atm_key_dct, ich_srt_mlf_bbn_keys)
    assert set(ich_srt_bbn_keys) == set(backbone_keys(rgr))
    atm_exp_hyd_keys_dct = atom_explicit_hydrogen_keys(rgr)
    last_num = len(ich_srt_bbn_keys)
    atm_ich_num_dct = {}
    for num, bbn_key in enumerate(ich_srt_bbn_keys):
        atm_ich_num_dct[bbn_key] = num
        for exp_hyd_key in atm_exp_hyd_keys_dct[bbn_key]:
            atm_ich_num_dct[exp_hyd_key] = last_num
            last_num += 1
    assert set(atm_ich_num_dct.keys()) == set(atm_keys)

    return ich, atm_ich_num_dct
