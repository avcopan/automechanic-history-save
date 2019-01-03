""" molecular graph -> InChI string conversion
"""
from .res import atom_bond_valences
from .res import atom_radical_valences
from .res import bond_orders
from ._shared import atom_keys
from ._shared import bond_keys
from ._shared import atom_symbols
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
    assert len(atm_keys) == len(nums)
    atm_ich_num_dct = dict(zip(atm_keys, nums))
    return ich, atm_ich_num_dct
