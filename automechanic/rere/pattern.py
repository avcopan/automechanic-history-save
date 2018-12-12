""" re pattern generators
"""
from re import escape as re_escape


def maybe(pattern):
    """ a pattern that may or may not be present

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return r'(?:{:s})?'.format(pattern)


def not_followed_by(pattern):
    """ a pattern is not present

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return r'(?!{:s})'.format(pattern)


def escape(pattern):
    """ escape special characters in pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return re_escape(pattern)


def zero_or_more(pattern, greedy=True):
    """ zero or more repeats of a pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param greedy: match as much as possible?
    :type greedy: bool

    :rtype: str
    """
    return (r'(?:{:s})*'.format(pattern) if greedy else
            r'(?:{:s})*?'.format(pattern))


def one_or_more(pattern, greedy=True):
    """ one or more repeats of a pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param greedy: match as much as possible?
    :type greedy: bool

    :rtype: str
    """
    return (r'(?:{:s})+'.format(pattern) if greedy else
            r'(?:{:s})+?'.format(pattern))


def capturing(pattern):
    """ generate a capturing pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return r'({:s})'.format(pattern)


def named_capturing(pattern, name):
    """ generate a capturing pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param name: a name for the capture
    :type name: str

    :rtype: str
    """
    return r'(?P<{:s}>{:s})'.format(name, pattern)


def one_of_these(patterns):
    """ any one of a series of patterns

    :param patterns: a series of `re` patterns
    :type patterns: list of strings

    :rtype: str
    """
    return r'(?:{:s})'.format('|'.join(patterns))
