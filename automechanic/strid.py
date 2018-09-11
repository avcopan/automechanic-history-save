""" string-ID-based functions
"""
from .ipybel.smiles import canonical as canonical_smiles
from .ipybel.smiles import number_of_atoms as number_of_atoms_from_smiles
from .ipybel.smiles import formula as formula_from_smiles
from .ipybel.smiles import geometry as geometry_from_smiles
from .ipybel.smiles import xyz_string as xyz_string_from_smiles
from .parse import DIGIT
from .parse import one_or_more
from .parse import named_capture
from .parse import group_dictionary


def canonical(sid):
    """ canonical SMILES string in a species ID
    """
    mult = multiplicity(sid)
    smi = smiles(sid)
    can_smi = canonical_smiles(smi)
    can_sid = '{:s}_m{:d}'.format(can_smi, mult)
    return can_sid


def smiles(sid):
    """ SMILES string from a species ID
    """
    smi = sid.split('_')[0]
    return smi


def multiplicity(sid):
    """ multiplicity from a species ID
    """
    mult_pattern = '_m' + named_capture(one_or_more(DIGIT), name='mult')
    gdct = group_dictionary(mult_pattern, sid)
    mult = int(gdct['mult'])
    return mult


def spin_count(sid):
    """ 2 * S = multiplicity - 1
    """
    return multiplicity(sid) - 1


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


def canonical_complex_identifier(cid):
    """ canonical SMILES for a set of species IDs sid1.sid2...
    """
    sids = cid.split('.')
    can_cid = '.'.join(sorted(map(canonical, sids)))
    return can_cid


def canonical_reaction_identifier(rid):
    """ canonical SMILES for a reaction identifier
    """
    rct_cid, prd_cid = rid.split('>>')
    can_rct_cid = canonical_complex_identifier(rct_cid)
    can_prd_cid = canonical_complex_identifier(prd_cid)
    can_rid = '>>'.join(sorted((can_rct_cid, can_prd_cid)))
    return can_rid


def split_reaction_identifier(rid):
    """ SMIRKS-style reaction ID from reactant and product species IDs
    """
    rct_str, prd_str = rid.split('>>')
    rct_sids = tuple(rct_str.split('.'))
    prd_sids = tuple(prd_str.split('.'))
    return (rct_sids, prd_sids)


def reaction_spin_counts(rid):
    """ multiplicities of reactions and products
    """
    rct_sids, prd_sids = split_reaction_identifier(rid)
    rct_mults = tuple(map(spin_count, rct_sids))
    prd_mults = tuple(map(spin_count, prd_sids))
    return (rct_mults, prd_mults)


def is_radical_radical(rid):
    """ determine if this is a radical-radical abstraction
    """
    ret = any((all(sct > 0 for sct in scts) and len(scts) > 1)
              for scts in reaction_spin_counts(rid))
    return ret


def is_spin_balanced(rid):
    """ determine if this reaction has equal total spin on both sides
    """
    rct_scts, prd_scts = reaction_spin_counts(rid)
    return sum(rct_scts) == sum(prd_scts)


if __name__ == '__main__':
    print canonical_reaction_identifier(
        'O[O]_m2.C=C[C](C(C=C)C)C_m2>>[O][O]_m3.C=CC(C(C=C)C)C_m1')
    print canonical_reaction_identifier(
        'C=CC(C(C=C)C)C_m1.[O][O]_m3>>C=C[C](C(C=C)C)C_m2.O[O]_m2')
