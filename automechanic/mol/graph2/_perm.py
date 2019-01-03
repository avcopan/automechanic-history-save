"""Permutation signature."""
from sympy.combinatorics.permutations import _af_parity


def parity(pmt, ref):
    """parity of a permutation

    :param pmt: permutation, as a sequence of elements
    :type pmt: tuple
    :param ref: identity permutation
    :type ref: tuple

    :rtype: bool (False = even, True = odd)
    """
    idxs = tuple(tuple(ref).index(e) for e in pmt)
    return bool(_af_parity(idxs))
