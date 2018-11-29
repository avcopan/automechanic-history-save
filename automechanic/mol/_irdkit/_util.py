""" common helpers
"""
import rdkit.Chem as _rd_chem


def inchi_from_rd_mol(rdm, options='', strict=True):
    """ InChI string from an rdkit molecule object
    """
    ich = _rd_chem.inchi.MolToInchi(
        rdm, options=options, treatWarningAsError=strict)
    return ich
