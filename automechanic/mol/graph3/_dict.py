""" dictionary helpers
"""


def by_key(dct, keys, fill_val=None):
    """ dictionary on a set of keys, filling missing entries
    """
    return dict(zip(keys, values_by_key(dct, keys, fill_val=fill_val)))


def values_by_key(dct, keys, fill_val=None):
    """ return dictionary values for specific keys, filling missing entries
    """
    return tuple(dct[key] if key in dct else fill_val for key in keys)


def transform_values(dct, func):
    """ apply a function to each value
    """
    return dict(zip(dct.keys(), map(func, dct.values())))
