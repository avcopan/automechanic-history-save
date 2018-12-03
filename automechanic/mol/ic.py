""" functions operating on InChI strings
"""
import itertools
from .ge import inchi as _inchi_from_geometry
from ._irdkit import from_inchi as _rdm_from_inchi
from ._irdkit import to_inchi as _rdm_to_inchi
from ._irdkit import geometry as _rdm_to_geometry
from ._irdkit import inchi_to_inchi_key as _inchi_to_inchi_key
from ..rere.pattern import escape as _escape
from ..rere.pattern import capturing as _capturing
from ..rere.pattern import zero_or_more as _zero_or_more
from ..rere.pattern import one_or_more as _one_or_more
from ..rere.pattern import one_of_these as _one_of_these
from ..rere.pattern_lib import UPPERCASE_LETTER as _UPPERCASE_LETTER
from ..rere.pattern_lib import NONWHITESPACE as _NONWHITESPACE
from ..rere.pattern_lib import STRING_START as _STRING_START
from ..rere.pattern_lib import STRING_END as _STRING_END
from ..rere.find import single_capture as _single_capture
from ..rere.find import ends_with as _ends_with
from ..rere.find import replace as _replace
from ..rere.find import split as _split


def inchi(ich, force_stereo=False):
    """ recompute InChI string
    """
    _options = '-SUU' if force_stereo else ''
    rdm = _rdm_from_inchi(ich)
    ich = _rdm_to_inchi(rdm, options=_options, with_aux_info=False)
    return ich


def inchi_prefix(ich):
    """ InChI prefix
    """
    _start = _STRING_START
    _body = _escape('InChI=') + _one_or_more(_NONWHITESPACE, greedy=False)
    _end = _escape('/')
    _pattern = _start + _capturing(_body) + _end
    return _single_capture(_pattern, ich)


def inchi_formula(ich):
    """ InChI formula
    """
    _start = (_STRING_START + _escape('InChI=') +
              _one_or_more(_NONWHITESPACE, greedy=False) + _escape('/'))
    _body = _UPPERCASE_LETTER + _zero_or_more(_NONWHITESPACE, greedy=False)
    _end = _escape('/')
    _pattern = _start + _capturing(_body) + _end
    return _single_capture(_pattern, ich)


def inchi_sublayer(ich, key):
    """ a sublayer from the InChI string, by key
    """
    slash = _escape('/')
    _start = slash + key
    _body = _one_or_more(_NONWHITESPACE, greedy=False)
    _end = _one_of_these([slash, _STRING_END])
    _pattern = _start + _capturing(_body) + _end
    return _single_capture(_pattern, ich)


def inchi_base(ich):
    """ InChI string with only the main layer
    """
    ich_pfx = inchi_prefix(ich)
    ich_fml = inchi_formula(ich)
    ich_con = 'c' + inchi_sublayer(ich, key='c')
    ich_hyd = 'h' + inchi_sublayer(ich, key='h')
    ich_bas = '/'.join([ich_pfx, ich_fml, ich_con, ich_hyd])
    return ich_bas


def normalized_inchi(ich):
    """ standard InChI string with main layer and stereo b and t layers
    """
    ich_bas = inchi_base(ich)

    parts = [ich_bas]
    ich_b_sub = inchi_sublayer(ich, key='b')
    if ich_b_sub:
        ich_ez = 'b' + ich_b_sub
        parts.append(ich_ez)

    ich_t_sub = inchi_sublayer(ich, key='t')
    if ich_t_sub:
        ich_th = 't' + ich_t_sub
        parts.append(ich_th)

    ich_ret = '/'.join(parts)
    return ich_ret


def _known_stereo_elements(ich, key):
    assert key in ('b', 't')
    ich_b = inchi_sublayer(ich, key)
    eles = _split(',', ich_b) if ich_b else ()
    return tuple(ele for ele in eles if not _ends_with(ele, _escape('?')))


def known_double_bond_stereo_elements(ich):
    """ known double bond stereo elements
    """
    return _known_stereo_elements(ich, 'b')


def known_tetrahedral_stereo_elements(ich):
    """ known tetrahedral stereo elements
    """
    return _known_stereo_elements(ich, 't')


def _expand_unknown_stereo_element(ele):
    _pattern = _escape('?') + _STRING_END
    if not _ends_with(_pattern, ele):
        ret = (ele,)
    else:
        ele_m = _replace(_pattern, '-', ele)
        ele_p = _replace(_pattern, '+', ele)
        ret = (ele_m, ele_p)
    return ret


def _all_double_bond_stereo_elements(ich):
    ich_fs = inchi(ich, force_stereo=True)
    ich_fs_b = inchi_sublayer(ich_fs, 'b')
    return _split(',', ich_fs_b) if ich_fs_b else None


def _all_tetrahedral_stereo_elements(ich):
    _pattern = _escape('?') + _STRING_END
    ich_fs = inchi(ich, force_stereo=True)
    ich_fs_t = inchi_sublayer(ich_fs, 't')
    ret = None
    if ich_fs_t:
        eles = _split(',', ich_fs_t)
        ele0 = eles[0]
        ele0 = _replace(_pattern, '-', ele0)
        ret = (ele0,) + tuple(eles[1:])
    return ret


def _double_bond_stereo_layer_from_elements(eles):
    return 'b' + ','.join(eles)


def _tetrahedral_stereo_layer_from_elements(eles):
    return 't' + ','.join(eles)


def inchi_expand_unknown_stereo(ich):
    """ expand InChI string into its stereomers
    """
    ich_bas = inchi_base(ich)
    ich_b_eles = _all_double_bond_stereo_elements(ich)
    ich_t_eles = _all_tetrahedral_stereo_elements(ich)

    ich_lst = [ich_bas]
    if ich_b_eles:
        ich_b_eles_lst = itertools.product(
            *map(_expand_unknown_stereo_element, ich_b_eles))
        ich_b_lst = list(map(_double_bond_stereo_layer_from_elements,
                             ich_b_eles_lst))
        ich_lst = ['/'.join((ich_start, ich_end)) for ich_start, ich_end in
                   itertools.product(ich_lst, ich_b_lst)]

    if ich_t_eles:
        ich_t_eles_lst = itertools.product(
            *map(_expand_unknown_stereo_element, ich_t_eles))
        ich_t_lst = list(map(_tetrahedral_stereo_layer_from_elements,
                             ich_t_eles_lst))
        ich_lst = ['/'.join((ich_start, ich_end)) for ich_start, ich_end in
                   itertools.product(ich_lst, ich_t_lst)]

    return tuple(ich_lst)


def inchi_key(ich):
    """ computes InChIKey from an InChI string
    """
    return _inchi_to_inchi_key(ich)


def geometry(ich):
    """ cartesian geometry from an InChI string
    """
    rdm = _rdm_from_inchi(ich)
    geo = _rdm_to_geometry(rdm)
    ich_ = _inchi_from_geometry(geo)

    # make sure the geometry matches the original inchi
    ich1 = inchi_base(ich)
    ich2 = inchi_base(ich_)
    assert inchi_key(ich1) == inchi_key(ich2)
    ez1 = known_double_bond_stereo_elements(ich)
    ez2 = known_double_bond_stereo_elements(ich_)
    assert set(ez1) <= set(ez2)
    th1 = known_tetrahedral_stereo_elements(ich)
    th2 = known_tetrahedral_stereo_elements(ich_)
    assert set(th1) <= set(th2)

    return geo
