""" functions operating on SMILES strings
"""
from ._irdkit import from_smiles as _rdm_from_smiles
from ._irdkit import to_inchi as _rdm_to_inchi


def inchi(smi, force_stereo=False):
    """ InChI string from a SMILES string
    """
    rdm = _rdm_from_smiles(smi)
    _options = '-SUU' if force_stereo else ''
    ich = _rdm_to_inchi(rdm, options=_options, with_aux_info=False)
    return ich
