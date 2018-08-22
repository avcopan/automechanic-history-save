""" a library of specialized SMILES functions
"""
from .smiles import number_of_atoms
from .smiles import formula
from .smiles import split_smirks
from .smiles import make_smirks
from .formlib import abstraction_argsort as formula_abstraction_argsort


def match_hydrogen_abstraction_formula(smrk):
    """ see if SMIRKS matches H abstraction; return sorted SMIRKS if so
    """
    rsmis, psmis = split_smirks(smrk)
    rfmls = tuple(map(formula, rsmis))
    pfmls = tuple(map(formula, psmis))
    abstr_srt = formula_abstraction_argsort(rfmls, pfmls)
    ret_smrk = None
    if abstr_srt:
        r_posns, p_posns = abstr_srt
        q1h_smi, q2_smi = (rsmis[posn] for posn in r_posns)
        q1_smi, q2h_smi = (psmis[posn] for posn in p_posns)
        if number_of_atoms(q1h_smi) >= number_of_atoms(q2h_smi):
            ret_smrk = make_smirks((q1h_smi, q2_smi), (q1_smi, q2h_smi))
        else:
            ret_smrk = make_smirks((q2h_smi, q1_smi), (q2_smi, q1h_smi))
    return ret_smrk


if __name__ == '__main__':
    SMRK = '[CH2]C=C.[CH2]C=C>>C=C=C.CC=C'
    print(match_hydrogen_abstraction_formula(SMRK))
