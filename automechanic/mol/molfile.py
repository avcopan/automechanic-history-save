""" functions operating on MOLFile V3000 strings
"""
from ._irdkit import from_molfile as _rdm_from_molfile
from ._irdkit import to_inchi as _rdm_to_inchi


def inchi(mlf):
    """ InChI string from a MOLFile string
    """
    rdm = _rdm_from_molfile(mlf)
    ich = _rdm_to_inchi(rdm)
    return ich
