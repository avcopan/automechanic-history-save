""" tuple dictionary helpers

dictionary values must all be iterables of the same length
"""
import numpy
from ._dict import transform_values as _transform_values


def by_key_by_position(dct, keys, pos):
    """ position values, as a dictionary with these keys
    """
    assert set(keys) <= set(dct.keys())
    vals = tuple(map(dct.__getitem__, keys))
    assert pos <= numpy.shape(vals)[1]
    pos_vals = list(zip(*vals))[pos]
    return dict(zip(keys, pos_vals))


def set_by_key_by_position(dct, pos_dct, pos):
    """ set values by position and key
    """
    assert set(pos_dct.keys()) <= set(dct.keys())
    dct = _transform_values(dct, func=list)
    for key, pos_val in pos_dct.items():
        dct[key][pos] = pos_val
    dct = _transform_values(dct, func=tuple)
    return dct
