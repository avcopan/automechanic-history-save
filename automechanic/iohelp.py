""" helpers for the io module
"""
from .geom import graph
from .geom import xyz_string
from .graph import atom_neighborhood_indices
from .strid import split_reaction_identifier
from .geomlib import abstraction_indices as geom_abstraction_indices
from .stridlib import match_hydrogen_abstraction_formula


def find_abstraction(rid, mgeo_dct):
    """
    """
    dxyzs = None

    if match_hydrogen_abstraction_formula(rid):
        (q1h_sid, q2_sid), (q1_sid, q2h_sid) = split_reaction_identifier(rid)

        q1h_mgeo, q2_mgeo, q1_mgeo, q2h_mgeo = map(
            mgeo_dct.__getitem__, [q1h_sid, q2_sid, q1_sid, q2h_sid])

        idxs = geom_abstraction_indices(q1h_mgeo, q2_mgeo, q1_mgeo, q2h_mgeo)

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
    RID = 'C=CC_m1.[OH]_m2>>[CH2]C(O)C_m1'

    R1_SID = 'C=CC_m1'
    R2_SID = '[OH]_m2'
    P_SID = '[CH2]C(O)C_m1'

    R1_MGEO = (('C', (0.999749614329, 0.0815463557276, -0.05313151319214)),
               ('C', (0.2941554428798, -0.1333084924908, 1.06457863301)),
               ('C', (-1.198418794122, -0.141000266906, 1.104594274387)),
               ('H', (0.50843512849, 0.2642722240302, -1.003700339629)),
               ('H', (2.084674708630, 0.0790932629613, -0.040370798372)),
               ('H', (0.817176204583, -0.313348598161, 2.00117588120)),
               ('H', (-1.559370502273, -1.112686909705, 1.454564921229)),
               ('H', (-1.559370055401, 0.631682276162, 1.789879903006)),
               ('H', (-1.63349683031, 0.0486453884978, 0.1180237523367)))
    R2_MGEO = (('O', (0.952632526586, -0.0627316796020, 0.0080078583889)),
               ('H', (1.892440257233, -0.0627316796020, 0.0080078583889)))
    P_MGEO = (('C', (0.977343848365, -0.0764188131994, 0.13889740679)),
              ('H', (0.4624939859419, -1.02959411007, 0.155086468595)),
              ('H', (0.4440163978138, 0.783646033730, 0.528720307800)),
              ('C', (2.455537814018, -0.01753155811344, -0.002585548076801)),
              ('O', (2.853647414718, -0.998191131877, -0.945748607211)),
              ('C', (2.93814855121, 1.339952372503, -0.48915762596)),
              ('H', (2.921436713659, -0.251407370808, 0.959743143408)),
              ('H', (2.55606620911, -1.86222435081, -0.612788411425)),
              ('H', (2.608102938152, 2.147148597200, 0.172659593754)),
              ('H', (4.03290647398, 1.356511004561, -0.532775456763)),
              ('H', (2.58295049411, 1.55110745037, -1.504509980310)))

    MGEO_DCT = {R1_SID: R1_MGEO,
                R2_SID: R2_MGEO,
                P_SID: P_MGEO}
