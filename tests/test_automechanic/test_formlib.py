""" test the automechanic.formlib module
"""
from itertools import permutations
from automechanic import formlib


def test__abstraction_argsort():
    """ test formlib.abstraction_argsort()
    """
    rct_fmls = ({'C': 1, 'H': 4}, {'O': 1, 'H': 1})
    prd_fmls = ({'C': 1, 'H': 3}, {'O': 1, 'H': 2})

    for rct_posns in permutations(range(2)):
        for prd_posns in permutations(range(2)):
            _rct_fmls = [rct_fmls[posn] for posn in rct_posns]
            _prd_fmls = [prd_fmls[posn] for posn in prd_posns]
            abstr_srt = formlib.abstraction_argsort(_rct_fmls, _prd_fmls)
            assert abstr_srt == (rct_posns, prd_posns)


if __name__ == '__main__':
    test__abstraction_argsort()
