""" functions operating on cartesian geometries (angstroms)
"""
import numpy
from ._irdkit import from_molfile as _rdm_from_molfile2
from ._irdkit import to_inchi as _rdm_to_inchi
from ..rere.pattern import escape as _escape
from ..rere.pattern import capturing as _capturing
from ..rere.pattern import zero_or_more as _zero_or_more
from ..rere.pattern_lib import UNSIGNED_INTEGER as _UNSIGNED_INTEGER
from ..rere.find import split as _split
from ..rere.find import first_capture as _first_capture


def inchi(geo):
    """ InChI string of a cartesian geometry
    """
    ich, _ = inchi_with_order(geo)
    return ich


def _parse_inchi_order_from_auxinfo(ich_aux):
    _comma = _escape(',')
    _pattern = _escape('/N:') + _capturing(
        _zero_or_more(_UNSIGNED_INTEGER + _comma) + _UNSIGNED_INTEGER)
    one_index_order_str = _first_capture(_pattern, ich_aux)
    one_index_order = tuple(map(int, _split(_comma, one_index_order_str)))
    order = tuple(numpy.subtract(one_index_order, 1))
    return order


def inchi_with_order(geo):
    """ InChI string of a cartesian geometry
    """
    mbl = molfile(geo)
    rdm = _rdm_from_molfile2(mbl)
    ich, ich_aux = _rdm_to_inchi(rdm, with_aux_info=True)
    ich_ord = _parse_inchi_order_from_auxinfo(ich_aux)
    return ich, ich_ord


def molfile(geo):
    """ MOLFile string of a cartesian geometry
    """
    # TODO: implement this
    raise NotImplementedError
