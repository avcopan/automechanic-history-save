""" some dictionary helpers
"""


def by_key(dct, keys, fill=None):
    """ return dictionary for a set of keys, filling missing entries
    """
    return dict(zip(keys, values_by_key(dct, keys, fill=fill)))


def values_by_key(dct, keys, fill=None):
    """ return dictionary values for specific keys, filling missing entries
    """
    return tuple(dct[key] if key in dct else fill for key in keys)


def filter_by_value(dct, func=lambda val: val):
    """ return dictionary for a set of values, defined by a func
    """
    return {key: val for key, val in dct.items() if func(val)}
