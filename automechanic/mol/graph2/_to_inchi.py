""" molecular graph -> InChI string conversion
"""
from .res import atom_bond_valences
from .res import atom_radical_valences
from .res import bond_orders
from ._shared import atom_keys
from ._shared import bond_keys
from ._shared import backbone_keys
from ._shared import atom_symbols
from ._shared import atom_explicit_hydrogen_keys
from ._shared import highspin_resonance_graph
from ._molfile import from_data as _mlf_from_data
from ._irdkit import from_molfile as _rdm_from_molfile
from ._irdkit import to_inchi as _rdm_to_inchi
from ._dict import values_by_key as _values_by_key
from .._inchi_aux import numbering as _ich_aux_numbering


def inchi(xgr):
    """ InChI string of a molecular graph (no stereo)
    """
    ich, _ = _inchi_with_atom_priorities(xgr)
    return ich


def atom_inchi_numbers(xgr):
    """ InChI numbering of backbone atoms
    """
    _, nums = _inchi_with_atom_priorities(xgr)
    return nums


def _inchi_with_atom_priorities(xgr):
    """ InChI string of this resonance graph
    """
    rgr = highspin_resonance_graph(xgr)
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
    atm_ich_num_dct = _inchi_numbering(nums, xgr)

    return ich, atm_ich_num_dct


def _inchi_numbering(nums, xgr):
    """ inchi numbering, augmented to account for explicit hydrogens

    explicit hydrogens are added at the end, grouped and sorted by their
    associated backbone atom
    """
    atm_exp_hyd_keys_dct = atom_explicit_hydrogen_keys(xgr)
    bbn_keys = backbone_keys(xgr)

    assert len(bbn_keys) == len(nums)
    atm_ich_num_dct = dict(zip(bbn_keys, nums))

    last_num = max(nums)
    for bbn_key in backbone_keys(xgr):
        for exp_hyd_key in atm_exp_hyd_keys_dct[bbn_key]:
            last_num += 1
            atm_ich_num_dct[exp_hyd_key] = last_num

    return atm_ich_num_dct
