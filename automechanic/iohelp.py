""" helpers for the io module
"""
from .geom import graph
from .geom import xyz_string
from .graph import atom_neighborhood_indices
from .strid import split_reaction_identifier
from .geomlib import find_abstraction as find_abstraction_from_geometries
from .stridlib import match_hydrogen_abstraction_formula


def find_abstraction(rid, mgeo_dct):
    """
    """
    dxyzs = None

    if match_hydrogen_abstraction_formula(rid):
        (q1h_sid, q2_sid), (q1_sid, q2h_sid) = split_reaction_identifier(rid)

        q1h_mgeo, q2_mgeo, q1_mgeo, q2h_mgeo = map(
            mgeo_dct.__getitem__, [q1h_sid, q2_sid, q1_sid, q2h_sid])

        idxs = find_abstraction_from_geometries(
            q1h_mgeo, q2_mgeo, q1_mgeo, q2h_mgeo)

        if idxs:
            q1h_idx, q2_idx, _, _ = idxs

            q1h_dxyz = abstractee_xyz_string(q1h_mgeo, q1h_idx)
            q2_dxyz = abstractor_xyz_string(q2_mgeo, q2_idx)
            q1_dxyz = xyz_string(q1_mgeo)
            q2h_dxyz = xyz_string(q2h_mgeo)

            dxyzs = q1h_dxyz, q2_dxyz, q1_dxyz, q2h_dxyz

    return dxyzs


def abstractor_xyz_string(q_mgeo, q_idx):
    """ .xyz string with a '4' in front X*...H-Y-Z
    """
    dxyz = _xyz_string_with_labels(q_mgeo, lbl_dct={q_idx: 4})
    return dxyz


def abstractee_xyz_string(qh_mgeo, qh_idx):
    """ .xyz string with '2', '1', '3' in front X...H*-Y*-Z*
    """
    if len(qh_mgeo) < 3:
        raise ValueError("Too few atoms to form abstractee .xyz string")
    qh_mgrph = graph(qh_mgeo)
    y_idx = next(iter(atom_neighborhood_indices(qh_mgrph, qh_idx)))
    z_idx = next(i for i in atom_neighborhood_indices(qh_mgrph, y_idx)
                 if i != qh_idx)
    dxyz = _xyz_string_with_labels(qh_mgeo,
                                   lbl_dct={qh_idx: 2, y_idx: 1, z_idx: 3})
    return dxyz


def _xyz_string_with_labels(mgeo, lbl_dct):
    """
    """
    natms = len(mgeo)
    dxyz = '{:d}\n\n'.format(natms)
    for idx, (asymb, xyz) in enumerate(mgeo):
        if idx in lbl_dct:
            dxyz += '{:s} '.format(repr(lbl_dct[idx]))
        dxyz += '{:s} {:s} {:s} {:s}\n'.format(asymb, *map(repr, xyz))
    return dxyz


if __name__ == '__main__':
    RID = 'CC=CCCC=CC_m1.O[O]_m2>>CC=C[CH]CC=CC_m2.OO_m1'

    Q1H_SID = 'CC=CCCC=CC_m1'
    Q2_SID = 'O[O]_m2'
    Q1_SID = 'CC=C[CH]CC=CC_m2'
    Q2H_SID = 'OO_m1'

    Q1H_MGEO = (('C', (0.94026, 0.06311, 0.03590)),
                ('C', (0.43390, 1.40773, 0.43771)),
                ('C', (-0.36222, 2.14515, -0.34927)),
                ('C', (-0.85534, 3.52342, -0.00798)),
                ('C', (-0.35711, 4.58091, -1.00204)),
                ('C', (-0.85022, 4.32541, -2.39876)),
                ('C', (-1.64634, 5.15644, -3.08615)),
                ('C', (-2.15270, 4.83847, -4.45303)),
                ('H', (0.60455, -0.21937, -0.96774)),
                ('H', (0.58750, -0.69409, 0.74082)),
                ('H', (2.03416, 0.06482, 0.03701)),
                ('H', (0.74621, 1.78901, 1.40600)),
                ('H', (-0.67271, 1.73926, -1.30976)),
                ('H', (-1.95215, 3.51055, 0.00839)),
                ('H', (-0.52756, 3.80739, 0.99942)),
                ('H', (-0.68488, 5.56884, -0.65638)),
                ('H', (0.73970, 4.59804, -1.01387)),
                ('H', (-0.53974, 3.39183, -2.86321)),
                ('H', (-1.95866, 6.09933, -2.64579)),
                ('H', (-3.24661, 4.83947, -4.45125)),
                ('H', (-1.79994, 5.58882, -5.16523)),
                ('H', (-1.81700, 3.85420, -4.79697)))
    Q2_MGEO = (('O', (1.09058, -0.08516, 0.00491)),
               ('O', (0.59345, 0.53149, 1.16172)),
               ('H', (2.06256, -0.08222, 0.01041)))
    Q1_MGEO = (('C', (1.03809, -0.04849, -0.12581)),
               ('C', (0.60859, 0.74811, 1.04286)),
               ('C', (-0.27191, 1.76442, 1.08468)),
               ('C', (-0.79235, 2.28643, 2.33648)),
               ('H', (-1.80506, 2.68585, 2.29059)),
               ('C', (-0.33989, 1.74579, 3.67133)),
               ('C', (-0.81457, 0.32236, 3.83828)),
               ('C', (0.00768, -0.74027, 3.88941)),
               ('C', (-0.45074, -2.15215, 3.75287)),
               ('H', (0.68361, 0.36638, -1.06832)),
               ('H', (0.65340, -1.07021, -0.01494)),
               ('H', (2.13022, -0.10284, -0.12813)),
               ('H', (1.04870, 0.37913, 1.95408)),
               ('H', (-0.70832, 2.15212, 0.16963)),
               ('H', (0.75339, 1.82220, 3.73876)),
               ('H', (-0.74047, 2.36360, 4.48138)),
               ('H', (-1.88788, 0.16658, 3.75942)),
               ('H', (1.08670, -0.59425, 3.90534)),
               ('H', (-0.03950, -2.58850, 2.83431)),
               ('H', (-0.09449, -2.73759, 4.60232)),
               ('H', (-1.54046, -2.22498, 3.69809)))
    Q2H_MGEO = (('O', (0.94095, 0.01468, 0.00913)),
                ('O', (0.56219, 0.91718, 1.08344)),
                ('H', (1.90630, 0.10786, 0.12005)),
                ('H', (-0.40316, 0.82400, 0.97252)))

    MGEO_DCT = {Q1H_SID: Q1H_MGEO,
                Q2_SID: Q2_MGEO,
                Q1_SID: Q1_MGEO,
                Q2H_SID: Q2H_MGEO}

    print find_abstraction(RID, MGEO_DCT)
