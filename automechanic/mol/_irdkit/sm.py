""" rdkit SMILES string interface
"""
import rdkit.Chem as _rd_chem
from ._util import inchi_from_rd_mol as _inchi_from_rd_mol


def inchi(smi, options='', strict=True):
    """ InChI string of a SMILES string
    """
    rdm = _rd_mol(smi)
    return _inchi_from_rd_mol(rdm, options=options, strict=strict)


def _rd_mol(smi):
    rdm = _rd_chem.MolFromSmiles(smi)
    return _rd_chem.AddHs(rdm)
