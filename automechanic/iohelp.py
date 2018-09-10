""" helpers for the io module
"""
from itertools import chain
from itertools import permutations
import pandas
from more_itertools import unique_everseen
from .strid import reaction_identifier
from .strid import split_reaction_identifier
from .strid import canonical_reaction_identifier
from .strid import formula as formula_from_strid
from .strid import number_of_atoms as number_of_atoms_from_strid
from .form import subtract as subtract_formulas
from .geom import graph
from .geom import xyz_string
from .geom import addition_indices
from .geom import abstraction_indices
from .geom import migration_indices
from .graph import atom_neighborhood_indices


def merge_reaction_dataframes(rxn_df1, rxn_df2):
    """ ugly function to merge reaction dataframes with 'reaction_id' columns
    """
    assert all('reaction_id' in df for df in (rxn_df1, rxn_df2))
    can_rid_iter1 = map(canonical_reaction_identifier,
                        rxn_df1['reaction_id'])
    can_rid_iter2 = map(canonical_reaction_identifier,
                        rxn_df2['reaction_id'])
    rxn_keys1 = list(frozenset(map(frozenset, split_reaction_identifier(rid)))
                     for rid in can_rid_iter1)
    rxn_keys2 = list(frozenset(map(frozenset, split_reaction_identifier(rid)))
                     for rid in can_rid_iter2)

    rxn_keys = set(rxn_keys1) | set(rxn_keys2)
    rows = []
    for rxn_key in rxn_keys:
        row = {}
        if rxn_key in rxn_keys1:
            idx1 = rxn_keys1.index(rxn_key)
            row1 = {key: val for key, val in rxn_df1.iloc[idx1].items()
                    if not pandas.isna(val)}
            row.update(row1)
        if rxn_key in rxn_keys2:
            idx2 = rxn_keys2.index(rxn_key)
            row2 = {key: val for key, val in rxn_df2.iloc[idx2].items()
                    if not pandas.isna(val)}
            row.update(row2)
        assert row
        rows.append(row)

    rxn_df = pandas.DataFrame(rows)
    cols = [col for col in unique_everseen(chain(rxn_df1.columns,
                                                 rxn_df2.columns))
            if col in rxn_df]
    return rxn_df[cols]


def addition_candidate(rid):
    """ identify addition candidates, based on some loose criteria
    """
    return bool(_sorted_addition_candidates(rid))


def abstraction_candidate(rid):
    """ identify abstraction candidates, based on some loose criteria
    """
    return bool(_sorted_abstraction_candidate(rid))


def migration_candidate(rid):
    """ identify migration candidate, based on some loose criteria
    """
    rct_sids, prd_sids = split_reaction_identifier(rid)
    return len(rct_sids) == len(prd_sids) == 1


def addition(rid, mgeo_dct):
    """ sorted reaction ID with abstraction indices (or None)
    """
    addtn = None
    can_rids = _sorted_addition_candidates(rid)
    for can_rid in can_rids:
        rct_mgeos, prd_mgeos = _reaction_geometries(can_rid, mgeo_dct)
        x_mgeo, y_mgeo = rct_mgeos
        xy_mgeo, = prd_mgeos
        idxs = addition_indices(x_mgeo, y_mgeo, xy_mgeo)
        if idxs:
            addtn = (can_rid, idxs)
    return addtn


def abstraction(rid, mgeo_dct):
    """ sorted reaction ID with abstraction indices (or None)
    """
    abstr = None
    can_rid = _sorted_abstraction_candidate(rid)
    if can_rid:
        rct_mgeos, prd_mgeos = _reaction_geometries(can_rid, mgeo_dct)
        q1h_mgeo, q2_mgeo = rct_mgeos
        q1_mgeo, q2h_mgeo = prd_mgeos
        idxs = abstraction_indices(q1h_mgeo, q2_mgeo, q1_mgeo, q2h_mgeo)
        if idxs:
            abstr = (can_rid, idxs)
    return abstr


def migration(rid, mgeo_dct):
    """ (trivially) sorted reaction ID with migration indices (or None)
    """
    mgrtn = None
    if migration_candidate(rid):
        rct_mgeos, prd_mgeos = _reaction_geometries(rid, mgeo_dct)
        r_mgeo, = rct_mgeos
        p_mgeo, = prd_mgeos
        idxs = migration_indices(r_mgeo, p_mgeo)
        if idxs:
            mgrtn = (rid, idxs)
    return mgrtn


def addition_xyz_strings(rid, idxs, mgeo_dct):
    """ TorsScan xyz strings for an addition reaction
    """
    dxyz_dct = None
    (x_sid, y_sid), (xy_sid,) = split_reaction_identifier(rid)
    (x_mgeo, y_mgeo), (xy_mgeo,) = _reaction_geometries(rid, mgeo_dct)
    x_idx, y_idx, _, _ = idxs
    xy_dxyz = xyz_string(xy_mgeo, {})
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
    q1h_idx, q2_idx, _, _ = idxs
    q1_dxyz = xyz_string(q1_mgeo, {})
    q2h_dxyz = xyz_string(q2h_mgeo, {})

    q1h_dxyz = _abstraction_xyz_string_q1h(q1h_mgeo, q1h_idx)
    q2_dxyz = _abstraction_xyz_string_q2(q2_mgeo, q2_idx)

    if q1h_dxyz:
        dxyz_dct = {q1_sid: q1_dxyz, q2h_sid: q2h_dxyz,
                    q1h_sid: q1h_dxyz, q2_sid: q2_dxyz}

    return dxyz_dct


def migration_xyz_strings(rid, idxs, mgeo_dct):
    """ TorsScan xyz strings for a migration reaction
    """
    dxyz_dct = None
    (r_sid,), (p_sid,) = split_reaction_identifier(rid)
    (r_mgeo,), (p_mgeo,) = _reaction_geometries(rid, mgeo_dct)

    r_idx_h, r_idx_a, _, _ = idxs

    p_dxyz = xyz_string(p_mgeo, {})
    r_dxyz = _migration_xyz_string(r_mgeo, r_idx_h, r_idx_a)

    if r_dxyz:
        dxyz_dct = {p_sid: p_dxyz, r_sid: r_dxyz}

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


def migration_input_string(rid, tmp_str, tmp_kevyal_dct):
    """ TorsScan input for migration reaction
    """
    (r_sid,), (p_sid,) = split_reaction_identifier(rid)
    sub_dct = {'r': r_sid, 'p': p_sid}
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
            (r1_sid, r2_sid, p1_sid, p2_sid)
            for r1_sid, r2_sid in permutations(rct_sids)
            for p1_sid, p2_sid in permutations(prd_sids)
            if _matches_abstraction_formula(r1_sid, r2_sid, p1_sid, p2_sid))
        sids = next(itr, None)
        if sids:
            q1h_sid, q2_sid, q1_sid, q2h_sid = sids
            if (number_of_atoms_from_strid(q1h_sid) < 3 and
                    number_of_atoms_from_strid(q2h_sid) >= 3):
                q1_sid, q2_sid = q2_sid, q1_sid
                q1h_sid, q2h_sid = q2h_sid, q1h_sid
            can_rid = reaction_identifier((q1h_sid, q2_sid), (q1_sid, q2h_sid))
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
    dxyz = xyz_string(mgeo, {y_idx: 4})
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
        dxyz = xyz_string(mgeo, {x_idx: 2, xn1_idx: 1, xn2_idx: 3})
    return dxyz


def _abstraction_xyz_string_q1h(mgeo, h_idx):
    # these functions happen to be the same
    return _addition_xyz_string_x(mgeo=mgeo, x_idx=h_idx)


def _abstraction_xyz_string_q2(mgeo, q_idx):
    # these functions happen to be the same
    return _addition_xyz_string_y(mgeo=mgeo, y_idx=q_idx)


def _migration_xyz_string(mgeo, h_idx, a_idx):
    dxyz = None
    mgrph = graph(mgeo)
    _2chainz = ((an1_idx, an2_idx)
                for an1_idx in atom_neighborhood_indices(mgrph, a_idx)
                for an2_idx in atom_neighborhood_indices(mgrph, an1_idx)
                if an2_idx != a_idx
                and an2_idx not in atom_neighborhood_indices(mgrph, h_idx))
    idxs = next(_2chainz, None)
    if idxs:
        an1_idx, an2_idx = idxs
        dxyz = xyz_string(mgeo, {h_idx: 1, a_idx: 2, an1_idx: 3, an2_idx: 4})
    return dxyz
