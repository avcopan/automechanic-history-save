""" generic regexes for the `re` module
"""
from re import escape


def maybe(pattern):
    """ a pattern that may or may not be present

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return r'(?:{:s})?'.format(pattern)


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


def repeat_range(pattern, nmin, nmax):
    """ `m` to `n` repeats of a pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param nmin: minimum number of repeats
    :type nmin: int
    :param nmax: maximum number of repeats
    :type nmax: int

    :rtype: str
    """
    return r'(?:{:s}){{{:d},{:d}}}'.format(pattern, nmin, nmax)


def capture(pattern):
    """ generate a capturing pattern

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return r'({:s})'.format(pattern)


def named_capture(pattern, name):
    """ generate a named capturing pattern

    :param pattern: an `re` pattern
    :type pattern: str
    :param name: a group name for the capture
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


STRING_START = r'\A'
STRING_END = r'\Z'
LINE_START = r'^'
LINE_END = r'$'
PERIOD = escape('.')
PLUS = escape('+')
MINUS = escape('-')
OPEN_PAREN = escape('(')
CLOSE_PAREN = escape(')')
LINE = r'^.*\n'
WHITESPACE = r'[ \t]'
WHITESPACES = one_or_more(WHITESPACE)
DIGIT = r'[0-9]'
INTEGER = one_or_more(DIGIT)
FLOAT = maybe(MINUS) + one_or_more(DIGIT) + PERIOD + one_or_more(DIGIT)
UPPERCASE_LETTER = r'[A-Z]'
LETTER = r'[a-zA-Z]'
