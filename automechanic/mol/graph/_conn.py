""" specific connnectivity graph functions
"""
from .to_inchi import with_numbering as _inchi_with_numbering


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
