""" helpers for the io module
"""
from itertools import permutations
from .pchemkin import split_reaction as split_chemkin_reaction
from .pchemkin import split_therm_data
from .strid import split_reaction_identifier
from .strid import formula as formula_from_strid
from .strid import reaction_identifier as reaction_identifier_from_sids
from .geom import graph
from .geom import xyz_string
from .geom import abstraction_indices
from .geom import addition_indices
from .geom import migration_indices
from .graph import indices as graph_indices
from .graph import atom_neighborhood_indices
from .form import subtract as subtract_formulas
from .therm import enthalpy as therm_enthalpy
from .timeout import timeout


def translate_chemkin_reaction(rxn_str, sid_dct):
    """ determine the reaction ID for a CHEMKIN reaction
    """
    rid = None

    spcs = sid_dct.keys()
    rcts, prds = split_chemkin_reaction(rxn_str)
    if set(rcts) | set(prds) < set(spcs):
        rct_sids = map(sid_dct.__getitem__, rcts)
        prd_sids = map(sid_dct.__getitem__, prds)
        rid = reaction_identifier_from_sids(rct_sids, prd_sids)

    return rid


def thermo_value_dictionary(thd_strs, sid_dct):
    """ a dictionary of thermo values
    """
    thv_dct = dict((sid_dct[spc], val) for spc, val in
                   map(translate_chemkin_thermo_data, thd_strs)
                   if spc in sid_dct)
    return thv_dct


def translate_chemkin_thermo_data(thd_str):
    """ determine the heat of formation from a CHEMKIN NASA polynomial
    """
    t_298 = 298.
    ret = None

    thd = split_therm_data(thd_str)

    if not thd:
        raise ValueError("Failed to parse thermo data:\n{:s}"
                         .format(str(thd_str)))

    spc, cfts_lo, cfts_hi, t_cross, t_lo, t_hi = thd

    assert t_lo < t_cross < t_hi
    cfts = cfts_lo if t_298 <= t_cross else cfts_hi
    h_298 = therm_enthalpy(t_298, cfts)
    ret = (spc, h_298)
    return ret


def abstraction_candidate(rid):
    """ identify abstraction candidates, based on some loose criteria
    """
    return bool(_abs_sorted_candidate(rid))


def addition_candidate(rid):
    """ identify addition candidates, based on some loose criteria
    """
    return bool(_add_sorted_candidates(rid))


def migration_candidate(rid):
    """ identify migration candidate, based on some loose criteria
    """
    return bool(_mig_sorted_candidate(rid))


@timeout(300)
def abstraction(rid, mgeo_dct, _thv_dct=None):
    """ species IDs with abstraction indices (or None)
    """
    abst = None
    abst_sids = _abs_sorted_candidate(rid)
    if abst_sids:
        mgeos = list(map(mgeo_dct.__getitem__, abst_sids))
        abst_idxs = abstraction_indices(*mgeos)
        if abst_idxs:
            if _thv_dct and _is_uphill(abst_sids[:2], abst_sids[2:], _thv_dct):
                abst_sids = tuple(reversed(abst_sids))
                abst_idxs = tuple(reversed(abst_idxs))
            if _abs_reverse_for_torsscan(abst_sids, abst_idxs, mgeo_dct):
                abst_sids = tuple(reversed(abst_sids))
                abst_idxs = tuple(reversed(abst_idxs))
            abst = (abst_sids, abst_idxs)
    return abst


@timeout(300)
def addition(rid, mgeo_dct, _thv_dct=None):
    """ species IDs with addition indices (or None)
    """
    addn = None
    addn_sids_lst = _add_sorted_candidates(rid)
    for addn_sids in addn_sids_lst:
        mgeos = list(map(mgeo_dct.__getitem__, addn_sids))
        addn_idxs = addition_indices(*mgeos)
        if addn_idxs:
            if _add_swap_for_torsscan(addn_sids, addn_idxs, mgeo_dct):
                addn_sids, addn_idxs = _add_switch_reactants(
                    addn_sids, addn_idxs)
            addn = (addn_sids, addn_idxs)
    return addn


@timeout(300)
def migration(rid, mgeo_dct, _thv_dct=None):
    """ species IDs with migration indices (or None)
    """
    mgrn = None
    mgrn_sids = _mig_sorted_candidate(rid)
    if mgrn_sids:
        mgeos = list(map(mgeo_dct.__getitem__, mgrn_sids))
        mgrn_idxs = migration_indices(*mgeos)
        if mgrn_idxs:
            if _thv_dct and _is_uphill(mgrn_sids[:1], mgrn_sids[1:], _thv_dct):
                mgrn_sids, mgrn_idxs = _mig_reverse(mgrn_sids, mgrn_idxs)
            mgrn = (mgrn_sids, mgrn_idxs)
    return mgrn


def abstraction_xyz_strings(sids, idxs, mgeo_dct):
    """ TorsScan xyz strings for an abstraction reaction
    """
    dxyz_dct = None
    q1h_sid, q2_sid, q1_sid, q2h_sid = sids
    q1h_mgeo, q2_mgeo, q1_mgeo, q2h_mgeo = map(mgeo_dct.__getitem__, sids)
    q1h_idx, q2_idx, _, _ = idxs
    q1_dxyz = xyz_string(q1_mgeo, {})
    q2h_dxyz = xyz_string(q2h_mgeo, {})

    q1h_dxyz = _abstraction_xyz_string_q1h(q1h_mgeo, q1h_idx)
    q2_dxyz = _abstraction_xyz_string_q2(q2_mgeo, q2_idx)

    if q1h_dxyz:
        dxyz_dct = {q1_sid: q1_dxyz, q2h_sid: q2h_dxyz,
                    q1h_sid: q1h_dxyz, q2_sid: q2_dxyz}

    return dxyz_dct


def abstraction_input_string(sids, tmp_str, tmp_kevyal_dct):
    """ TorsScan input for abstraction reaction
    """
    q1h_sid, q2_sid, q1_sid, q2h_sid = sids
    sub_dct = {'q1h': q1h_sid, 'q2': q2_sid, 'q1': q1_sid, 'q2h': q2h_sid}
    sub_dct.update(tmp_kevyal_dct)
    inp_str = tmp_str.format(**sub_dct)
    return inp_str


def addition_xyz_strings(sids, idxs, mgeo_dct):
    """ TorsScan xyz strings for an addition reaction
    """
    dxyz_dct = None
    x_sid, y_sid, xy_sid = sids
    x_mgeo, y_mgeo, xy_mgeo = map(mgeo_dct.__getitem__, sids)
    x_idx, y_idx, _, _ = idxs
    xy_dxyz = xyz_string(xy_mgeo, {})
    x_dxyz = _addition_xyz_string_x(x_mgeo, x_idx)
    y_dxyz = _addition_xyz_string_y(y_mgeo, y_idx)

    if x_dxyz:
        dxyz_dct = {xy_sid: xy_dxyz, x_sid: x_dxyz, y_sid: y_dxyz}

    return dxyz_dct


def addition_input_string(sids, tmp_str, tmp_kevyal_dct):
    """ TorsScan input for addition reaction
    """
    x_sid, y_sid, xy_sid = sids
    sub_dct = {'x': x_sid, 'y': y_sid, 'xy': xy_sid}
    sub_dct.update(tmp_kevyal_dct)
    inp_str = tmp_str.format(**sub_dct)
    return inp_str


def migration_xyz_strings(sids, idxs, mgeo_dct):
    """ TorsScan xyz strings for a migration reaction
    """
    dxyz_dct = None
    r_sid, p_sid = sids
    r_mgeo, p_mgeo = map(mgeo_dct.__getitem__, sids)

    r_idx_h, r_idx_a, _, _ = idxs

    p_dxyz = xyz_string(p_mgeo, {})
    r_dxyz = _migration_xyz_string(r_mgeo, r_idx_h, r_idx_a)

    if r_dxyz:
        dxyz_dct = {p_sid: p_dxyz, r_sid: r_dxyz}

    return dxyz_dct


def migration_input_string(sids, tmp_str, tmp_kevyal_dct):
    """ TorsScan input for migration reaction
    """
    r_sid, p_sid = sids
    sub_dct = {'r': r_sid, 'p': p_sid}
    sub_dct.update(tmp_kevyal_dct)
    inp_str = tmp_str.format(**sub_dct)
    return inp_str


# helpers
def _is_uphill(rct_sids, prd_sids, thv_dct):
    rct_thv = sum(map(thv_dct.__getitem__, rct_sids))
    prd_thv = sum(map(thv_dct.__getitem__, prd_sids))
    return rct_thv < prd_thv


def _abs_sorted_candidate(rid):
    sids = None
    rct_sids, prd_sids = split_reaction_identifier(rid)
    if len(rct_sids) == len(prd_sids) == 2:
        itr = (
            (r1_sid, r2_sid, p1_sid, p2_sid)
            for r1_sid, r2_sid in permutations(rct_sids)
            for p1_sid, p2_sid in permutations(prd_sids)
            if _abs_matches_formula(r1_sid, r2_sid, p1_sid, p2_sid))
        sids = next(itr, None)
    return sids


def _add_sorted_candidates(rid):
    sids_tup = ()
    prd_sids, rct_sids = sorted(split_reaction_identifier(rid), key=len)
    if len(rct_sids) == 2 and len(prd_sids) == 1:
        xy_sid, = prd_sids
        sids_tup = tuple((x_sid, y_sid, xy_sid) for x_sid, y_sid in
                         permutations(rct_sids))
    return sids_tup


def _mig_sorted_candidate(rid):
    sids = None
    rct_sids, prd_sids = split_reaction_identifier(rid)
    if len(rct_sids) == len(prd_sids) == 1:
        rct_sid, = rct_sids
        prd_sid, = prd_sids
        sids = (rct_sid, prd_sid)
    return sids


def _abs_matches_formula(r1_sid, r2_sid, p1_sid, p2_sid):
    r1_fml, r2_fml, p1_fml, p2_fml = map(
        formula_from_strid, (r1_sid, r2_sid, p1_sid, p2_sid))
    match = (subtract_formulas(p1_fml, r1_fml) == {'H': -1} and
             subtract_formulas(p2_fml, r2_fml) == {'H': +1})
    return match


def _abs_reverse_for_torsscan(abst_sids, abst_idxs, mgeo_dct):
    q1h_sid, _, _, q2h_sid = abst_sids
    q1h_idx, _, _, q2h_idx = abst_idxs
    q1h_mgeo = mgeo_dct[q1h_sid]
    q2h_mgeo = mgeo_dct[q2h_sid]
    ret = (not _abstraction_xyz_string_q1h(q1h_mgeo, q1h_idx) and
           _abstraction_xyz_string_q1h(q2h_mgeo, q2h_idx))
    return ret


def _add_swap_for_torsscan(addn_sids, addn_idxs, mgeo_dct):
    x_sid, y_sid, _ = addn_sids
    x_idx, y_idx, _, _ = addn_idxs
    x_mgeo = mgeo_dct[x_sid]
    y_mgeo = mgeo_dct[y_sid]
    ret = (not _addition_xyz_string_x(x_mgeo, x_idx) and
           _addition_xyz_string_x(y_mgeo, y_idx))
    return ret


def _mig_reverse(mgrn_sids, mgrn_idxs):
    assert len(mgrn_sids) == 2
    assert len(mgrn_idxs) == 4
    rev_mgrn_sids = tuple(mgrn_sids[1:] + mgrn_sids[:1])
    rev_mgrn_idxs = tuple(mgrn_idxs[2:] + mgrn_idxs[:2])
    return rev_mgrn_sids, rev_mgrn_idxs


def _add_switch_reactants(addn_sids, addn_idxs):
    assert len(addn_sids) == 3
    assert len(addn_idxs) == 4
    x_sid, y_sid, xy_sid = addn_sids
    x_idx, y_idx, xy_idx_x, xy_idx_y = addn_idxs
    ret_addn_sids = (y_sid, x_sid, xy_sid)
    ret_addn_idxs = (y_idx, x_idx, xy_idx_y, xy_idx_x)
    return (ret_addn_sids, ret_addn_idxs)


def _addition_xyz_string_x(mgeo, x_idx):
    dxyz = None
    mgrph = graph(mgeo)
    _2chainz = ((xn1_idx, xn2_idx)
                for xn1_idx in atom_neighborhood_indices(mgrph, x_idx)
                for xn2_idx in atom_neighborhood_indices(mgrph, xn1_idx)
                if xn2_idx != x_idx)
    idxs = next(_2chainz, None)

    # allow for broken chains, if not enough connected atoms
    if not idxs:
        _broken_chains = ((x1_idx, x2_idx)
                          for x1_idx in atom_neighborhood_indices(mgrph, x_idx)
                          for x2_idx in graph_indices(mgrph)
                          if x2_idx not in (x_idx, x1_idx))
        idxs = next(_broken_chains, None)

    if idxs:
        xn1_idx, xn2_idx = idxs
        dxyz = xyz_string(mgeo, {x_idx: 2, xn1_idx: 1, xn2_idx: 3})
    return dxyz


def _addition_xyz_string_y(mgeo, y_idx):
    dxyz = xyz_string(mgeo, {y_idx: 4})
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
