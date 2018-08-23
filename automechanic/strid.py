""" string-ID-based functions
"""
from .ipybel.smiles import number_of_atoms as number_of_atoms_from_smiles
from .ipybel.smiles import formula as formula_from_smiles
from .ipybel.smiles import geometry as geometry_from_smiles
from .ipybel.smiles import xyz_string as xyz_string_from_smiles


def smiles(sid):
    """ SMILES string from a species ID
    """
    smi = sid.split('_')[0]
    return smi


def formula(sid):
    """ molecular formula
    """
    smi = smiles(sid)
    fml = formula_from_smiles(smi)
    return fml


def geometry(sid):
    """ molecular geometry
    """
    smi = smiles(sid)
    mgeo = geometry_from_smiles(smi)
    return mgeo


def xyz_string(sid):
    """ .xyz string
    """
    smi = smiles(sid)
    dxyz = xyz_string_from_smiles(smi)
    return dxyz


def number_of_atoms(sid):
    """ number of atoms
    """
    smi = smiles(sid)
    fml = number_of_atoms_from_smiles(smi)
    return fml


def reaction_identifier(rct_sids, prd_sids):
    """ SMIRKS-style reaction ID from reactant and product species IDs
    """
    rct_str = '.'.join(rct_sids)
    prd_str = '.'.join(prd_sids)
    rid = rct_str + '>>' + prd_str
    return rid


def split_reaction_identifier(rid):
    """ SMIRKS-style reaction ID from reactant and product species IDs
    """
    rct_str, prd_str = rid.split('>>')
    rct_sids = tuple(rct_str.split('.'))
    prd_sids = tuple(prd_str.split('.'))
    return (rct_sids, prd_sids)
