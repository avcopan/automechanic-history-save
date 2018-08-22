""" a library of specialized formula functions
"""
from itertools import permutations
from .form import subtract


def abstraction_argsort(rct_fmls, prd_fmls):
    """ find posible hydrogen abstraction, by reagent position
    """
    abstr_posn = None
    if len(rct_fmls) == len(prd_fmls) == 2:
        for rct_posns in permutations(range(2)):
            for prd_posns in permutations(range(2)):
                _rct_fmls = [rct_fmls[posn] for posn in rct_posns]
                _prd_fmls = [prd_fmls[posn] for posn in prd_posns]
                if matches_abstraction(_rct_fmls, _prd_fmls):
                    abstr_posn = rct_posns, prd_posns
    return abstr_posn


def matches_abstraction(rct_fmls, prd_fmls):
    """ reaction formula satisfies XHn + YHm -> XHn-1 + YHm+1 (exactly) ?
    """
    ret = False
    if len(rct_fmls) == len(prd_fmls) == 2:
        rct1_fml, rct2_fml = rct_fmls
        prd1_fml, prd2_fml = prd_fmls
        ret = (subtract(prd1_fml, rct1_fml) == {'H': -1} and
               subtract(prd2_fml, rct2_fml) == {'H': +1})
    return ret
