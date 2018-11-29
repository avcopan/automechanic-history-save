""" rdkit InChI string interface
"""
import rdkit.Chem as _rd_chem
import rdkit.Chem.AllChem as _rd_all_chem
from ._util import inchi_from_rd_mol as _inchi_from_rd_mol


def inchi(ich, options='', strict=True):
    """ recompute InChI string from an InChI string
    """
    rdm = _rd_mol(ich, strict=strict)
    return _inchi_from_rd_mol(rdm, options=options, strict=strict)


def inchi_key(ich, strict=True):
    """ computes InChIKey from an InChI string
    """
    ich = inchi(ich, strict=True) if strict else ich
    return _rd_chem.inchi.InchiToInchiKey(ich)


def atom_symbols_and_coordinates(ich, strict=True):
    """ cartesian geometry (angstroms) from an InChI string
    """
    rdm = _rd_mol(ich, strict=strict)
    _rd_all_chem.EmbedMolecule(rdm)
    asbs = tuple(rda.GetSymbol() for rda in rdm.GetAtoms())
    xyzs = tuple(map(tuple, rdm.GetConformer(0).GetPositions()))
    return (asbs, xyzs)


def _rd_mol(ich, strict=True):
    rdm = _rd_chem.inchi.MolFromInchi(ich, treatWarningAsError=strict)
    return _rd_chem.AddHs(rdm)
