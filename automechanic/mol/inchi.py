""" functions operating on InChI strings
"""
import itertools
import numpy
from .geom import inchi as _inchi_from_geometry
from .graph2.conn import (make_hydrogens_implicit as
                          _graph_conn_make_hydrogens_implicit)
from ._irdkit import from_inchi as _rdm_from_inchi
from ._irdkit import to_inchi as _rdm_to_inchi
from ._irdkit import to_smiles as _rdm_to_smiles
from ._irdkit import to_molfile as _rdm_to_molfile
from ._irdkit import inchi_to_inchi_key as _inchi_to_inchi_key
from ._irdkit import geometry as _rdm_to_geometry
from ._irdkit import connectivity_graph as _rdm_to_connectivity_graph
from ._ipybel import from_inchi as _pbm_from_inchi
from ._ipybel import geometry as _pbm_to_geometry
from ..rere.pattern import escape as _escape
from ..rere.pattern import capturing as _capturing
from ..rere.pattern import named_capturing as _named_capturing
from ..rere.pattern import zero_or_more as _zero_or_more
from ..rere.pattern import one_or_more as _one_or_more
from ..rere.pattern import one_of_these as _one_of_these
from ..rere.pattern import not_followed_by as _not_followed_by
# from ..rere.pattern_lib import ANY_CHAR as _ANY_CHAR
from ..rere.pattern_lib import LOWERCASE_LETTER as _LOWERCASE_LETTER
from ..rere.pattern_lib import UNSIGNED_INTEGER as _UNSIGNED_INTEGER
from ..rere.pattern_lib import NONWHITESPACE as _NONWHITESPACE
from ..rere.pattern_lib import STRING_START as _STRING_START
from ..rere.pattern_lib import STRING_END as _STRING_END
from ..rere.find import first_capture as _first_capture
from ..rere.find import first_named_capture as _first_named_capture
from ..rere.find import ends_with as _ends_with
from ..rere.find import replace as _replace
from ..rere.find import split as _split

_INCHI_SUBLAYER_END = _one_of_these([_escape('/'), _STRING_END])


class PARSE():
    """ InChI format specifications """

    class PREFIX():
        """ _ """
        _START = _STRING_START + _escape('InChI=')
        _VERSION = _one_or_more(_NONWHITESPACE, greedy=False)
        _END = _INCHI_SUBLAYER_END

        ALL_KEY = 'all'
        VERSION_KEY = 'version'
        PATTERN = _named_capturing(
            _START + _named_capturing(_VERSION, name=VERSION_KEY),
            name=ALL_KEY) + _END

    class FORMULA():
        """ _ """
        _START = _escape('/') + _not_followed_by(_LOWERCASE_LETTER)
        _FORMULA = _one_or_more(_NONWHITESPACE, greedy=False)
        _END = _INCHI_SUBLAYER_END

        ALL_KEY = 'all'
        PATTERN = _START + _named_capturing(_FORMULA, name=ALL_KEY) + _END


def smiles(ich):
    """ SMILES string from an InChI string
    """
    rdm = _rdm_from_inchi(ich)
    smi = _rdm_to_smiles(rdm)
    return smi


def molfile(ich):
    """ MOLFile string from an InChI string
    """
    rdm = _rdm_from_inchi(ich)
    mfl = _rdm_to_molfile(rdm)
    return mfl


def recalculate(ich, force_stereo=False):
    """ recalculate InChI string
    """
    _options = '-SUU' if force_stereo else ''
    rdm = _rdm_from_inchi(ich)
    ich = _rdm_to_inchi(rdm, options=_options, with_aux_info=False)
    return ich


def is_closed(ich):
    """ regenerating the InChI string yields the same thing
    """
    return recalculate(ich) == ich


def prefix(ich):
    """ InChI prefix
    """
    cap_dct = _first_named_capture(PARSE.PREFIX.PATTERN, ich)
    assert cap_dct
    pfx = cap_dct[PARSE.PREFIX.ALL_KEY]
    return pfx


def version(ich):
    """ InChI version
    """
    cap_dct = _first_named_capture(PARSE.PREFIX.PATTERN, ich)
    assert cap_dct
    ver = cap_dct[PARSE.PREFIX.VERSION_KEY]
    return ver


def formula_layer(ich):
    """ InChI formula
    """
    cap_dct = _first_named_capture(PARSE.FORMULA.PATTERN, ich)
    assert cap_dct
    fml = cap_dct[PARSE.FORMULA.ALL_KEY]
    return fml


def sublayer(ich, key):
    """ a sublayer from the InChI string, by key
    """
    slash = _escape('/')
    _start = slash + key
    _body = _one_or_more(_NONWHITESPACE, greedy=False)
    _end = _one_of_these([slash, _STRING_END])
    _pattern = _start + _capturing(_body) + _end
    return _first_capture(_pattern, ich)


def with_sublayers(ich, keys):
    """ InChI string with only the main layer
    """
    ich_pfx = prefix(ich)
    ich_fml = formula_layer(ich)

    parts = [ich_pfx, ich_fml]

    for key in keys:
        sub = sublayer(ich, key=key)
        if sub:
            parts.append('{key:s}{sub:s}'.format(key=key, sub=sub))

    ich_bas = '/'.join(parts)
    return ich_bas


def core_parent(ich):
    """ get the InChI string of the core parent structure
    """
    ich_cp = with_sublayers(ich, ('c', 'h'))
    assert is_closed(ich_cp)
    return ich_cp


def inchi_key(ich):
    """ computes InChIKey from an InChI string
    """
    return _inchi_to_inchi_key(ich)


def connectivity_graph(ich):
    """ connectivity graph from an InChI string
    """
    rdm = _rdm_from_inchi(ich)

    # make sure the InChI string was valid and that the graph will be
    # inchi-sorted
    ich_, ich_aux = _rdm_to_inchi(rdm, with_aux_info=True)
    ich_ord = _parse_inchi_order_from_auxinfo(ich_aux)
    assert list(ich_ord) == sorted(ich_ord)
    assert core_parent(ich) == core_parent(ich_)

    cgr = _rdm_to_connectivity_graph(rdm)
    cgr = _graph_conn_make_hydrogens_implicit(cgr)
    return cgr


def _parse_inchi_order_from_auxinfo(ich_aux):
    _comma = _escape(',')
    _pattern = _escape('/N:') + _capturing(
        _zero_or_more(_UNSIGNED_INTEGER + _comma) + _UNSIGNED_INTEGER)
    one_index_order_str = _first_capture(_pattern, ich_aux)
    one_index_order = tuple(map(int, _split(_comma, one_index_order_str)))
    order = tuple(numpy.subtract(one_index_order, 1))
    return order


# def stereo_graph(ich):
#     """ stereo graph from an InChI string
#     """
#     cgr = connectivity_graph(ich)
#     assert not has_unknown_stereo_elements(ich)
#     print(cgr)
#     _atom_stereo_values(ich)
#     _bond_stereo_values(ich)
#
#
# def _bond_stereo_values(ich):
#     eles = bond_stereo_elements(ich)
#     print(eles)


def geometry(ich):
    """ cartesian geometry from an InChI string
    """
    try:
        rdm = _rdm_from_inchi(ich)
        geo = _rdm_to_geometry(rdm)
        assert _inchi_matches_geometry(ich, geo)
    except (AssertionError, RuntimeError):
        pbm = _pbm_from_inchi(ich)
        geo = _pbm_to_geometry(pbm)
        assert _inchi_matches_geometry(ich, geo)
    return geo


def _inchi_matches_geometry(ich, geo):
    ich_from_geo = _inchi_from_geometry(geo)
    ich1 = core_parent(ich)
    ich2 = core_parent(ich_from_geo)
    assert ich1 == ich2
    ez1 = known_bond_stereo_elements(ich)
    ez2 = known_bond_stereo_elements(ich_from_geo)
    th1 = known_atom_stereo_elements(ich)
    th2 = known_atom_stereo_elements(ich_from_geo)
    ret = set(ez1) <= set(ez2) and set(th1) <= set(th2)
    return ret


def has_unknown_stereo_elements(ich):
    """ does this InChI have unknown stereo elements
    """
    return bool(unknown_bond_stereo_elements(ich) or
                unknown_atom_stereo_elements(ich))


def bond_stereo_elements(ich):
    """ double bond stereo elements
    """
    ich = recalculate(ich, force_stereo=True)
    sub = sublayer(ich, 'b')
    eles = _split(',', sub) if sub else ()
    return tuple(eles)


def atom_stereo_elements(ich):
    """ atom stereo elements
    """
    ich = recalculate(ich, force_stereo=True)
    sub = sublayer(ich, 't')
    eles = _split(',', sub) if sub else ()
    return tuple(eles)


def unknown_bond_stereo_elements(ich):
    """ known double bond stereo elements
    """
    return tuple(ele for ele in bond_stereo_elements(ich)
                 if _ends_with(_escape('?'), ele))


def unknown_atom_stereo_elements(ich):
    """ known atom stereo elements
    """
    return tuple(ele for ele in atom_stereo_elements(ich)
                 if _ends_with(_escape('?'), ele))


def known_bond_stereo_elements(ich):
    """ known double bond stereo elements
    """
    return tuple(ele for ele in bond_stereo_elements(ich)
                 if not _ends_with(_escape('?'), ele))


def known_atom_stereo_elements(ich):
    """ known atom stereo elements
    """
    return tuple(ele for ele in atom_stereo_elements(ich)
                 if not _ends_with(_escape('?'), ele))


def compatible_stereoisomers(ich):
    """ expand InChI string into its stereomers
    """
    def _bond_stereo_layer_from_elements(eles):
        return 'b' + ','.join(eles)

    def _atom_stereo_layer_from_elements(eles):
        return 't' + ','.join(eles)

    def _expand_unknown_stereo_element(ele):
        _pattern = _escape('?') + _STRING_END
        if not _ends_with(_pattern, ele):
            ret = (ele,)
        else:
            ele_m = _replace(_pattern, '-', ele)
            ele_p = _replace(_pattern, '+', ele)
            ret = (ele_m, ele_p)
        return ret

    ich_bas = core_parent(ich)
    ich_b_eles = list(bond_stereo_elements(ich))
    ich_t_eles = list(atom_stereo_elements(ich))

    # the first atom stereo center is always assigned to '-'
    if ich_t_eles:
        _pattern = _escape('?') + _STRING_END
        ich_t_eles[0] = _replace(_pattern, '-', ich_t_eles[0])

    ich_lst = [ich_bas]
    if ich_b_eles:
        ich_b_eles_lst = itertools.product(
            *map(_expand_unknown_stereo_element, ich_b_eles))
        ich_b_lst = list(map(_bond_stereo_layer_from_elements,
                             ich_b_eles_lst))
        ich_lst = ['/'.join((ich_start, ich_end)) for ich_start, ich_end in
                   itertools.product(ich_lst, ich_b_lst)]

    if ich_t_eles:
        ich_t_eles_lst = itertools.product(
            *map(_expand_unknown_stereo_element, ich_t_eles))
        ich_t_lst = list(map(_atom_stereo_layer_from_elements,
                             ich_t_eles_lst))
        ich_lst = ['/'.join((ich_start, ich_end)) for ich_start, ich_end in
                   itertools.product(ich_lst, ich_t_lst)]

    # recalculating restores the /m and /s sublayers
    ich_lst = tuple(map(recalculate, ich_lst))
    return ich_lst
