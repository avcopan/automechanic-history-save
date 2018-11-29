""" functions operating on SMILES strings
"""
from ._irdkit.sm import inchi as _inchi


def inchi(smi, force_stereo=False, strict=True):
    """ InChI string from a SMILES string
    """
    _options = '-SUU' if force_stereo else ''
    ich = _inchi(smi, options=_options, strict=strict)
    return ich
