""" rdkit interface
"""
from rdkit import RDLogger
import rdkit.Chem as _rd_chem

_LOGGER = RDLogger.logger()
_LOGGER.setLevel(RDLogger.ERROR)


def from_molfile(mfl):
    """ rdkit molecule object from a mol block string
    """
    rdm = _rd_chem.rdmolfiles.MolFromMolBlock(mfl, removeHs=False)
    assert rdm is not None
    return rdm


def to_inchi(rdm, with_aux_info=False):
    """ InChI string from an rdkit molecule object
    """
    if with_aux_info:
        ret = _rd_chem.inchi.MolToInchiAndAuxInfo(rdm)
    else:
        ret = _rd_chem.inchi.MolToInchi(rdm)
    return ret
