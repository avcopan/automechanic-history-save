""" functions operating on InChI-AuxInfo strings
"""
import numpy
from ..rere.pattern import escape as _escape
from ..rere.pattern import one_or_more as _one_or_more
from ..rere.pattern import one_of_these as _one_of_these
from ..rere.pattern import named_capturing as _named_capturing
from ..rere.pattern_lib import UNSIGNED_INTEGER as _UNSIGNED_INTEGER
from ..rere.pattern_lib import NONWHITESPACE as _NONWHITESPACE
from ..rere.pattern_lib import STRING_END as _STRING_END
from ..rere.find import first_named_capture as _first_named_capture
from ..rere.find import all_captures as _all_captures

_NONWHITESPACES_NONGREEDY = _one_or_more(_NONWHITESPACE, greedy=False)
_INCHI_SUBLAYER_END = _one_of_these([_escape('/'), _STRING_END])


class PARSE():
    """ InChI-AuxInfo format specifications """
    class NUMBERING():
        """ _ """
        LAYER_KEY = 'all'
        CONTENT_KEY = 'content'

        _START = _escape('/')
        _LAYER = (_escape('N:') +
                  _named_capturing(_NONWHITESPACES_NONGREEDY,
                                   name=CONTENT_KEY))
        _END = _INCHI_SUBLAYER_END

        PATTERN = _START + _named_capturing(_LAYER, name=LAYER_KEY) + _END

        class NUMBER():
            """ _ """
            PATTERN = _UNSIGNED_INTEGER


def numbering(ich_aux):
    """ zero-indexed numbering
    """
    cap_dct = _first_named_capture(PARSE.NUMBERING.PATTERN, ich_aux)
    lyr = cap_dct[PARSE.NUMBERING.CONTENT_KEY]
    one_index_nums = tuple(
        map(int, _all_captures(PARSE.NUMBERING.NUMBER.PATTERN, lyr)))
    nums = tuple(numpy.subtract(one_index_nums, 1))
    return nums
