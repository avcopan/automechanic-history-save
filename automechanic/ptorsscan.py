""" parsers for TorsScan files
"""
import re
from .parse import maybe
from .parse import one_or_more
from .parse import capture
from .parse import SIGN
from .parse import INTEGER
from .parse import FLOAT
from .parse import SPACE


def arrhenius(plog_str):
    """ get arrhenius constants from a rate.plog file
    """
    ret = None

    spaces = one_or_more(SPACE, greedy=False)
    exp = FLOAT + 'E' + maybe(SIGN) + INTEGER
    arrh_pattern = spaces.join(['REACS=Capture', capture(exp), capture(FLOAT),
                                capture(FLOAT)])
    finds = re.findall(arrh_pattern, plog_str)
    if finds:
        arrh_a, arrh_b, arrh_e = map(float, finds[0])
        ret = (arrh_a, arrh_b, arrh_e)

    return ret
