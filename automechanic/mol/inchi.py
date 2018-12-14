""" functions operating on InChI strings
"""
import itertools
from .geom import inchi as _inchi_from_geometry
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
from ..rere.pattern import one_or_more as _one_or_more
from ..rere.pattern import one_of_these as _one_of_these
from ..rere.pattern import not_followed_by as _not_followed_by
from ..rere.pattern_lib import LOWERCASE_LETTER as _LOWERCASE_LETTER
from ..rere.pattern_lib import NONWHITESPACE as _NONWHITESPACE
from ..rere.pattern_lib import STRING_START as _STRING_START
from ..rere.pattern_lib import STRING_END as _STRING_END
from ..rere.find import first_capture as _first_capture
from ..rere.find import ends_with as _ends_with
from ..rere.find import replace as _replace
from ..rere.find import split as _split


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
    _start = _STRING_START
    _body = _escape('InChI=') + _one_or_more(_NONWHITESPACE, greedy=False)
    _end = _escape('/')
    _pattern = _start + _capturing(_body) + _end
    return _first_capture(_pattern, ich)


def formula_layer(ich):
    """ InChI formula
    """
    _start = (_escape('/') + _not_followed_by(_LOWERCASE_LETTER))
    _body = _one_or_more(_NONWHITESPACE, greedy=False)
    _end = _one_of_these([_escape('/'), _STRING_END])
    _pattern = _start + _capturing(_body) + _end
    return _first_capture(_pattern, ich)


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


def connectivity_graph(ich):
    """ connectivity graph from an InChI string
    """
    rdm = _rdm_from_inchi(ich)
    cgr = _rdm_to_connectivity_graph(rdm)
    return cgr


def _inchi_matches_geometry(ich, geo):
    ich_from_geo = _inchi_from_geometry(geo)
    ich1 = core_parent(ich)
    ich2 = core_parent(ich_from_geo)
    assert ich1 == ich2
    ez1 = known_double_bond_stereo_elements(ich)
    ez2 = known_double_bond_stereo_elements(ich_from_geo)
    th1 = known_tetrahedral_stereo_elements(ich)
    th2 = known_tetrahedral_stereo_elements(ich_from_geo)
    ret = set(ez1) <= set(ez2) and set(th1) <= set(th2)
    return ret


def has_unknown_stereo_elements(ich):
    """ does this InChI have unknown stereo elements
    """
    return bool(unknown_double_bond_stereo_elements(ich) or
                unknown_tetrahedral_stereo_elements(ich))


def double_bond_stereo_elements(ich):
    """ double bond stereo elements
    """
    ich = recalculate(ich, force_stereo=True)
    sub = sublayer(ich, 'b')
    eles = _split(',', sub) if sub else ()
    return tuple(eles)


def tetrahedral_stereo_elements(ich):
    """ tetrahedral stereo elements
    """
    ich = recalculate(ich, force_stereo=True)
    sub = sublayer(ich, 't')
    eles = _split(',', sub) if sub else ()
    return tuple(eles)


def unknown_double_bond_stereo_elements(ich):
    """ known double bond stereo elements
    """
    return tuple(ele for ele in double_bond_stereo_elements(ich)
                 if _ends_with(_escape('?'), ele))


def unknown_tetrahedral_stereo_elements(ich):
    """ known tetrahedral stereo elements
    """
    return tuple(ele for ele in tetrahedral_stereo_elements(ich)
                 if _ends_with(_escape('?'), ele))


def known_double_bond_stereo_elements(ich):
    """ known double bond stereo elements
    """
    return tuple(ele for ele in double_bond_stereo_elements(ich)
                 if not _ends_with(_escape('?'), ele))


def known_tetrahedral_stereo_elements(ich):
    """ known tetrahedral stereo elements
    """
    return tuple(ele for ele in tetrahedral_stereo_elements(ich)
                 if not _ends_with(_escape('?'), ele))


def compatible_stereoisomers(ich):
    """ expand InChI string into its stereomers
    """
    def _double_bond_stereo_layer_from_elements(eles):
        return 'b' + ','.join(eles)

    def _tetrahedral_stereo_layer_from_elements(eles):
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
    ich_b_eles = list(double_bond_stereo_elements(ich))
    ich_t_eles = list(tetrahedral_stereo_elements(ich))

    # the first tetrahedral stereo center is always assigned to '-'
    if ich_t_eles:
        _pattern = _escape('?') + _STRING_END
        ich_t_eles[0] = _replace(_pattern, '-', ich_t_eles[0])

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

    # recalculating restores the /m and /s sublayers
    ich_lst = tuple(map(recalculate, ich_lst))
    return ich_lst
