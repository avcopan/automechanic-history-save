""" stereo graph -> InChI string conversion
"""
from .molfile import from_data as _mlf_from_data
from .intco import atom_stereo_coordinates
from ._rdkit import from_molfile as _rdm_from_molfile
from ._rdkit import to_inchi_with_aux_info as _rdm_to_inchi_with_aux_info
from .._base import atom_keys
from .._base import bond_keys
from .._base import atom_symbols
from .._base import atom_bond_valences
from .._base import atom_radical_valences
from .._base import bond_orders
from .._base import explicit_stereo_sites
from .._base import is_chiral
from .._dict import values_by_key as _values_by_key


def inchi(xgr):
    """ InChI string of this stereo graph
    """
    xgr = explicit_stereo_sites(xgr)
    atm_keys = atom_keys(xgr)
    bnd_keys = bond_keys(xgr)
    atm_syms = _values_by_key(atom_symbols(xgr), atm_keys)
    atm_bnd_vlcs = _values_by_key(atom_bond_valences(xgr), atm_keys)
    atm_rad_vlcs = _values_by_key(atom_radical_valences(xgr), atm_keys)
    is_chi = is_chiral(xgr)
    atm_xyzs = _values_by_key(atom_stereo_coordinates(xgr), atm_keys)
    bnd_ords = _values_by_key(bond_orders(xgr), bnd_keys)
    mlf, _ = _mlf_from_data(atm_keys, bnd_keys, atm_syms, atm_bnd_vlcs,
                            atm_rad_vlcs, bnd_ords,
                            is_chi=is_chi, atm_xyzs=atm_xyzs)
    rdm = _rdm_from_molfile(mlf)
    ich, _ = _rdm_to_inchi_with_aux_info(rdm)
    return ich
