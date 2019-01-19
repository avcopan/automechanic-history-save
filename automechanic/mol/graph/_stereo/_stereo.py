""" stereo graph -> InChI string conversion
"""
from ._intco import atom_stereo_coordinates
from .._base import is_chiral
from .._base import explicit_stereo_sites
from ..to_inchi import with_numbering as _inchi_with_numbering


def inchi(sgr):
    """ InChI string of this stereo graph
    """
    sgr = explicit_stereo_sites(sgr)
    is_chi = is_chiral(sgr)
    atm_xyz_dct = atom_stereo_coordinates(sgr)
    ich, _ = _inchi_with_numbering(sgr, is_chi=is_chi, atm_xyz_dct=atm_xyz_dct)
    return ich
