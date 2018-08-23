""" a library of specialized species ID functions
"""
from .strid import formula
from .strid import number_of_atoms
from .strid import reaction_identifier
from .strid import split_reaction_identifier
from .formlib import abstraction_argsort as abstraction_argsort_from_formulas


def match_hydrogen_abstraction_formula(rid):
    """ see if reaction ID matches H abstraction; return sorted reaction ID if so
    """
    ret_rid = None

    rct_sids, prd_sids = split_reaction_identifier(rid)
    rct_fmls = tuple(map(formula, rct_sids))
    prd_fmls = tuple(map(formula, prd_sids))

    abstr_srt = abstraction_argsort_from_formulas(rct_fmls, prd_fmls)
    if abstr_srt:
        rct_posns, prd_posns = abstr_srt
        q1h_sid, q2_sid = (rct_sids[posn] for posn in rct_posns)
        q1_sid, q2h_sid = (prd_sids[posn] for posn in prd_posns)
        if number_of_atoms(q1h_sid) >= number_of_atoms(q2h_sid):
            ret_rid = reaction_identifier((q1h_sid, q2_sid), (q1_sid, q2h_sid))
        else:
            ret_rid = reaction_identifier((q2h_sid, q1_sid), (q2_sid, q1h_sid))

    return ret_rid


if __name__ == '__main__':
    print match_hydrogen_abstraction_formula(
        '[H][H]_m1.[O]_m3>>[H]_m2.[OH]_m2')
