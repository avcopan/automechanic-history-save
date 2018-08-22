""" formula-based functions
"""


def subtract(fml1, fml2):
    """ the difference between two molecular formulas
    """
    atms = set(fml1) | set(fml2)
    diffs = {}
    for atm in atms:
        cnt1 = fml1[atm] if atm in fml1 else 0
        cnt2 = fml2[atm] if atm in fml2 else 0
        if cnt1 != cnt2:
            diffs[atm] = cnt1 - cnt2
    return diffs
