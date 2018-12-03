""" re patterns
"""
from .pattern import maybe
from .pattern import escape
from .pattern import one_of_these
from .pattern import one_or_more
from .pattern import zero_or_more

STRING_START = r'\A'
STRING_END = r'\Z'

LINE_START = r'^'
LINE_END = r'$'

ANY_CHAR = r'[\s\S]'
NONNEWLINE = r'[^\n]'
NEWLINE = r'\n'
NONWHITESPACE = r'\S'
WHITESPACE = r'\s'
NONNEWLINE_WHITESPACE = r'[ \t]'

UPPERCASE_LETTER = '[A-Z]'
LOWERCASE_LETTER = '[a-z]'

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

NUMERICAL_EXPONENT = one_of_these(['E', 'e']) + INTEGER
EXPONENTIAL_INTEGER = INTEGER + NUMERICAL_EXPONENT
EXPONENTIAL_FLOAT = FLOAT + NUMERICAL_EXPONENT
