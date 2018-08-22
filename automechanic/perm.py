""" permutation helpers
"""
from functools import reduce


def compose(*pmts):
    """ compose several permutations
    """
    return reduce(product, pmts)


def product(pmt1, pmt2):
    """ the product of two permutations
    """
    nelem = len(pmt1)
    assert sorted(pmt1) == sorted(pmt2) == list(range(nelem))
    prod = tuple(pmt2[i] for i in pmt1)
    return prod


def inverse(pmt):
    """ get the inverse of a permutation
    """
    nelem = len(pmt)
    assert sorted(pmt) == list(range(nelem))
    pmt = tuple(pmt)
    inv = tuple(pmt.index(i) for i in range(nelem))
    return inv
