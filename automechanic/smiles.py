""" SMILES-based functions
"""
from .ipybel.smiles import number_of_atoms
from .ipybel.smiles import formula
from .ipybel.smiles import geometry
from .ipybel.smiles import xyz_string


def reverse_smirks(smrk):
    """ reverse a SMIRKS reaction
    """
    rct_str, prd_str = smrk.split('>>')
    ret_smrk = '>>'.join((prd_str, rct_str))
    return ret_smrk


def split_smirks(smrk):
    """ split a SMIRKS string into reactant and product SMILES
    """
    rct_str, prd_str = smrk.split('>>')
    rsmis = tuple(rct_str.split('.'))
    psmis = tuple(prd_str.split('.'))
    return (rsmis, psmis)


def make_smirks(rsmis, psmis):
    """ join reactant and product SMILES to form a SMIRKS string
    """
    rct_str = '.'.join(rsmis)
    prd_str = '.'.join(psmis)
    smrk = rct_str + '>>' + prd_str
    return smrk


__all__ = ['number_of_atoms', 'formula', 'geometry', 'xyz_string']
