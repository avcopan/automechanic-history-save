""" functions operating on InChIKeys
"""
from ..rere.pattern_lib import STRING_START as _STRING_START
from ..rere.pattern_lib import STRING_END as _STRING_END
from ..rere.pattern_lib import UPPERCASE_LETTER as _UPPERCASE_LETTER
from ..rere.pattern import escape as _escape
from ..rere.pattern import named_capturing as _named_capturing
from ..rere.find import has_match as _has_match
from ..rere.find import first_named_capture as _first_named_capture

_HASH1 = _UPPERCASE_LETTER * 14
_HASH2 = _UPPERCASE_LETTER * 8
_SVP = _UPPERCASE_LETTER * 2 + _escape('-') + _UPPERCASE_LETTER

HASH1_KEY = 'hash1'
HASH2_KEY = 'hash2'
SVP_KEY = 'svp'
ICK_PATTERN = (
    _STRING_START +
    _named_capturing(_HASH1, name=HASH1_KEY) + _escape('-') +
    _named_capturing(_HASH2, name=HASH2_KEY) +
    _named_capturing(_SVP, name=SVP_KEY) + _STRING_END)
SVP_STANDARD_NEUTRAL = 'SA-N'


def is_valid(ick):
    """ is this a valid InChIKey?
    """
    assert isinstance(ick, (str, bytes, bytearray))
    return _has_match(ICK_PATTERN, ick)


def is_standard_neutral(ick):
    """ is this a standard, netural InChIKey?
    """
    assert is_valid(ick)
    cap_dct = _first_named_capture(ICK_PATTERN, ick)
    svp = cap_dct[SVP_KEY]
    return svp == SVP_STANDARD_NEUTRAL


def first_hash(ick):
    """ the first hash block, indicating connectivity
    """
    assert is_valid(ick)
    cap_dct = _first_named_capture(ICK_PATTERN, ick)
    hash1 = cap_dct[HASH1_KEY]
    return hash1


def second_hash(ick):
    """ the second hash block, indicating connectivity
    """
    assert is_valid(ick)
    cap_dct = _first_named_capture(ICK_PATTERN, ick)
    hash2 = cap_dct[HASH2_KEY]
    return hash2
