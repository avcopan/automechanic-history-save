""" helpers for the io module
"""
from itertools import permutations
from .strid import reaction_identifier
from .strid import split_reaction_identifier
from .geom import graph
from .graph import atom_neighborhood_indices
from .form import subtract as subtract_formulas
from .strid import formula as formula_from_strid
from .geomlib import addition_indices as addition_indices_from_geometries
from .geomlib import abstraction_indices as abstraction_indices_from_geometries


def addition_candidate(rid):
    """ identify addition candidates, based on some loose criteria
    """
    return bool(_sorted_addition_candidates(rid))


def abstraction_candidate(rid):
    """ identify abstraction candidates, based on some loose criteria
    """
    return bool(_sorted_abstraction_candidate(rid))


def addition_indices(rid, mgeo_dct):
    """ sorted reaction ID with abstraction indices (or None)
    """
    addtn = None
    can_rids = _sorted_addition_candidates(rid)
    for can_rid in can_rids:
        rct_mgeos, prd_mgeos = _reaction_geometries(can_rid, mgeo_dct)
        x_mgeo, y_mgeo = rct_mgeos
        xy_mgeo, = prd_mgeos
        idxs = addition_indices_from_geometries(x_mgeo, y_mgeo, xy_mgeo)
        if idxs:
            addtn = (can_rid, idxs)
    return addtn


def abstraction_indices(rid, mgeo_dct):
    """ sorted reaction ID with abstraction indices (or None)
    """
    abstr = None
    can_rid = _sorted_abstraction_candidate(rid)
    if can_rid:
        rct_mgeos, prd_mgeos = _reaction_geometries(can_rid, mgeo_dct)
        q1h_mgeo, q2_mgeo = rct_mgeos
        q1_mgeo, q2h_mgeo = prd_mgeos
        idxs = abstraction_indices_from_geometries(q1h_mgeo, q2_mgeo,
                                                   q1_mgeo, q2h_mgeo)
        if idxs:
            abstr = (can_rid, idxs)
    return abstr


def addition_xyz_strings(rid, idxs, mgeo_dct):
    """ TorsScan xyz strings for an addition reaction
    """
    dxyz_dct = None
    (x_sid, y_sid), (xy_sid,) = split_reaction_identifier(rid)
    (x_mgeo, y_mgeo), (xy_mgeo,) = _reaction_geometries(rid, mgeo_dct)
    x_idx, y_idx = idxs
    xy_dxyz = _labeled_xyz_string(xy_mgeo, {})
    x_dxyz = _addition_xyz_string_x(x_mgeo, x_idx)
    y_dxyz = _addition_xyz_string_y(y_mgeo, y_idx)

    if x_dxyz:
        dxyz_dct = {xy_sid: xy_dxyz, x_sid: x_dxyz, y_sid: y_dxyz}

    return dxyz_dct


def abstraction_xyz_strings(rid, idxs, mgeo_dct):
    """ TorsScan xyz strings for an abstraction reaction
    """
    dxyz_dct = None
    (q1h_sid, q2_sid), (q1_sid, q2h_sid) = split_reaction_identifier(rid)
    (q1h_mgeo, q2_mgeo), (q1_mgeo, q2h_mgeo) = _reaction_geometries(rid,
                                                                    mgeo_dct)
    q1h_idx, q2_idx = idxs
    q1_dxyz = _labeled_xyz_string(q1_mgeo, {})
    q2h_dxyz = _labeled_xyz_string(q2h_mgeo, {})

    q1h_dxyz = _abstraction_xyz_string_q1h(q1h_mgeo, q1h_idx)
    q2_dxyz = _abstraction_xyz_string_q2(q2_mgeo, q2_idx)

    if q1h_dxyz:
        dxyz_dct = {q1_sid: q1_dxyz, q2h_sid: q2h_dxyz,
                    q1h_sid: q1h_dxyz, q2_sid: q2_dxyz}

    return dxyz_dct


def addition_input_string(rid, tmp_str, tmp_kevyal_dct):
    """ TorsScan input for addition reaction
    """
    (x_sid, y_sid), (xy_sid,) = split_reaction_identifier(rid)
    sub_dct = {'x': x_sid, 'y': y_sid, 'xy': xy_sid}
    sub_dct.update(tmp_kevyal_dct)
    inp_str = tmp_str.format(**sub_dct)
    return inp_str


def abstraction_input_string(rid, tmp_str, tmp_kevyal_dct):
    """ TorsScan input for abstraction reaction
    """
    (q1h_sid, q2_sid), (q1_sid, q2h_sid) = split_reaction_identifier(rid)
    sub_dct = {'q1h': q1h_sid, 'q2': q2_sid, 'q1': q1_sid, 'q2h': q2h_sid}
    sub_dct.update(tmp_kevyal_dct)
    inp_str = tmp_str.format(**sub_dct)
    return inp_str


def _sorted_addition_candidates(rid):
    can_rids = ()
    prd_sids, rct_sids = sorted(split_reaction_identifier(rid), key=len)
    if len(rct_sids) == 2 and len(prd_sids) == 1:
        xy_sid, = prd_sids
        can_rids = tuple(reaction_identifier((x_sid, y_sid), (xy_sid,))
                         for x_sid, y_sid in permutations(rct_sids))
    return can_rids


def _sorted_abstraction_candidate(rid):
    can_rid = None
    rct_sids, prd_sids = split_reaction_identifier(rid)
    if len(rct_sids) == len(prd_sids) == 2:
        itr = (
            reaction_identifier((r1_sid, r2_sid), (p1_sid, p2_sid))
            for r1_sid, r2_sid in permutations(rct_sids)
            for p1_sid, p2_sid in permutations(prd_sids)
            if _matches_abstraction_formula(r1_sid, r2_sid, p1_sid, p2_sid))
        can_rid = next(itr, None)
    return can_rid


def _matches_abstraction_formula(r1_sid, r2_sid, p1_sid, p2_sid):
    r1_fml, r2_fml, p1_fml, p2_fml = map(
        formula_from_strid, (r1_sid, r2_sid, p1_sid, p2_sid))
    match = (subtract_formulas(p1_fml, r1_fml) == {'H': -1} and
             subtract_formulas(p2_fml, r2_fml) == {'H': +1})
    return match


def _reaction_geometries(rid, mgeo_dct):
    rct_sids, prd_sids = split_reaction_identifier(rid)
    rct_mgeos = tuple(map(mgeo_dct.__getitem__, rct_sids))
    prd_mgeos = tuple(map(mgeo_dct.__getitem__, prd_sids))
    return rct_mgeos, prd_mgeos


def _addition_xyz_string_y(mgeo, y_idx):
    dxyz = _labeled_xyz_string(mgeo, {y_idx: 4})
    return dxyz


def _addition_xyz_string_x(mgeo, x_idx):
    dxyz = None
    mgrph = graph(mgeo)
    _2chainz = ((xn1_idx, xn2_idx)
                for xn1_idx in atom_neighborhood_indices(mgrph, x_idx)
                for xn2_idx in atom_neighborhood_indices(mgrph, xn1_idx)
                if xn2_idx != x_idx)
    idxs = next(_2chainz, None)
    if idxs:
        xn1_idx, xn2_idx = idxs
        dxyz = _labeled_xyz_string(mgeo, {x_idx: 2, xn1_idx: 1, xn2_idx: 3})
    return dxyz


def _abstraction_xyz_string_q1h(mgeo, h_idx):
    # these functions happen to be the same
    return _addition_xyz_string_x(mgeo=mgeo, x_idx=h_idx)


def _abstraction_xyz_string_q2(mgeo, q_idx):
    # these functions happen to be the same
    return _addition_xyz_string_y(mgeo=mgeo, y_idx=q_idx)


def _labeled_xyz_string(mgeo, lbl_dct):
    natms = len(mgeo)
    dxyz = '{:d}\n\n'.format(natms)
    for idx, (asymb, xyz) in enumerate(mgeo):
        if idx in lbl_dct:
            dxyz += '{:s} '.format(repr(lbl_dct[idx]))
        dxyz += '{:s} {:s} {:s} {:s}\n'.format(asymb, *map(repr, xyz))
    return dxyz
