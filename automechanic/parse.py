""" generic regexes for the `re` module
"""
from re import search
from re import finditer
from re import escape
from re import MULTILINE


# pattern makers
def not_followed_by(pattern):
    """ a pattern to exclude after a match

    :param pattern: an `re` pattern
    :type pattern: str

    :rtype: str
    """
    return r'(?!{:s})'.format(pattern)


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


# finders
def group_dictionary(pattern, string, flags=MULTILINE):
    """ return the group dictionary for an `re` search

    :param pattern: an `re` pattern with a named capture
    :type pattern: str
    :param string: a string
    :type string: str

    :rtype: str
    """
    gdct = None
    match = search(pattern, string, flags=flags)
    if match and match.groupdict():
        gdct = match.groupdict()
    return gdct


def group_lists(pattern, string, flags=MULTILINE):
    """ return group lists for all matches in an `re` finditer

    :param pattern: an `re` pattern with a named capture
    :type pattern: str
    :param string: a string
    :type string: str

    :rtype: str
    """
    glsts = [match.groups()
             for match in finditer(pattern, string, flags=flags)]
    return glsts


ANY_CHAR = r'[\s\S]'

NEWLINE = r'\n'
NON_NEWLINE = r'[^\n]'

SPACE = r'[ \t]'
SPACE_OR_NEWLINE = r'.'

UPPERCASE_LETTER = r'[A-Z]'
LETTER = r'[a-zA-Z]'

DIGIT = r'[0-9]'

PLUS = escape('+')
MINUS = escape('-')
SIGN = one_of_these([PLUS, MINUS])
UNSIGNED_INTEGER = one_or_more(DIGIT)
INTEGER = maybe(SIGN) + UNSIGNED_INTEGER

PERIOD = escape('.')
UNSIGNED_FLOAT = one_of_these(
    [zero_or_more(DIGIT) + PERIOD + one_or_more(DIGIT),
     one_or_more(DIGIT) + PERIOD + zero_or_more(DIGIT)])
FLOAT = maybe(SIGN) + UNSIGNED_FLOAT

STRING_START = r'\A'
STRING_END = r'\Z'
LINE_START = r'^'
LINE_END = r'$'
LINE = r'^.*\n'
