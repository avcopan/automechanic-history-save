""" rdkit mol block string interface
"""
import rdkit.Chem as _rd_chem


def inchi_with_auxinfo(mbl, options='', strict=True):
    """ InChI and InChI AuxInfo of a mol block string
    """
    rdm = _rdkit_mol(mbl)
    ich, ich_aux = _rd_chem.inchi.MolToInchiAndAuxInfo(
        rdm, options=options, treatWarningAsError=strict)
    return ich, ich_aux


def _rdkit_mol(mbl):
    return _rd_chem.rdmolfiles.MolFromMolBlock(mbl, removeHs=False)
