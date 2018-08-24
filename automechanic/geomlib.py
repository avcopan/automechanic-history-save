""" a library of specialized geom functions
"""
from itertools import product
from .geom import graph
from .geom import resonance_graphs
from .graph import atom_neighborhood_indices
from .graphlib import abstraction_indices as graph_abstraction_indices


def find_hydrogen_abstraction(r1h_mgeo, r2_mgeo, p1_mgeo, p2h_mgeo):
    """ find hydrogen abstraction for a pre-sorted reaction
    """
    idxs = None
    idxs1 = abstraction_indices(r1h_mgeo, p1_mgeo)
    idxs2 = abstraction_indices(p2h_mgeo, r2_mgeo)
    if idxs1 and idxs2:
        q1h_idx, q1_idx = idxs1
        q2h_idx, q2_idx = idxs2
        idxs = (q1h_idx, q2_idx, q1_idx, q2h_idx)
    return idxs


def abstraction_indices(qh_mgeo, q_mgeo):
    """ identify abstraction indices
    """
    idxs = None
    for qh_mgrph, q_mgrph in product(resonance_graphs(qh_mgeo),
                                     resonance_graphs(q_mgeo)):
        idxs = graph_abstraction_indices(qh_mgrph, q_mgrph)
        if idxs:
            break
    return idxs


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
    R1_MGEO = (('C', (0.94334996769, 0.0249179323630, -0.072794170880)),
               ('C', (0.43157430032, 1.40605117254, 0.36054528551)),
               ('C', (0.94334426836, 2.47190181665, -0.61887557282)),
               ('C', (-1.10375170209, 1.40604399655, 0.360551264139)),
               ('C', (0.94334877136, 1.72132671665, 1.77331403235)),
               ('H', (0.59366410943, -0.22884214957, -1.07983357206)),
               ('H', (0.59366781521, -0.75865490865, 0.60876436412)),
               ('H', (2.0386721840, -0.0062895482321, -0.082590942095)),
               ('H', (0.59365654557, 2.2734396369, -1.638252844)),
               ('H', (0.59365667194, 3.47090119292, -0.33512060245)),
               ('H', (2.03866510740, 2.4959915838, -0.64101216551)),
               ('H', (-1.4996964326, 0.65368274831, 1.05190269129)),
               ('H', (-1.49970101471, 2.3809584599, 0.66643845017)),
               ('H', (-1.4997001466, 1.18349687226, -0.6366947887)),
               ('H', (0.59366715174, 0.97609107155, 2.4965931307)),
               ('H', (2.03867112901, 1.72845675313, 1.8052372017)),
               ('H', (0.59366187178, 2.70336673358, 2.11112837508)))
    R2_MGEO = (('O', (0.92407436880, 0.0046449910356, 0.0311054435523)),
               ('O', (2.2060743689, 0.0046449910356, 0.0311054435523)))
    P1_MGEO = (('C', (1.04944543638, -0.07824486600, 0.036328012776)),
               ('C', (0.54215536886, 1.36214719030, 0.197067778203)),
               ('C', (1.08427456498, 2.22746354885, -0.91116608306)),
               ('H', (0.60559993711, 3.1759054439, -1.1353951895)),
               ('H', (2.11652774748, 2.11149339267, -1.2273793231)),
               ('C', (-0.99296520488, 1.3605915670, 0.160665400256)),
               ('C', (1.03174989061, 1.94200860225, 1.52919010147)),
               ('H', (0.69984189193, -0.52390691182, -0.9022819812)),
               ('H', (0.69295239563, -0.71379114437, 0.85615333921)),
               ('H', (2.14461430573, -0.123175741162, 0.036694067204)),
               ('H', (-1.40455935135, 0.76386009081, 0.983846103)),
               ('H', (-1.40019738061, 2.37406826899, 0.252496722990)),
               ('H', (-1.3712297983, 0.93512045333, -0.77619897532)),
               ('H', (0.6677139938, 1.35191290718, 2.37796336260)),
               ('H', (2.12702101473, 1.95532435005, 1.57778346608)),
               ('H', (0.6831102991, 2.97252820577, 1.66568733456)))
    P2_MGEO = (('O', (1.09603557697, -0.086284073525, 0.039206622023)),
               ('O', (0.59891144765, 1.22442998853, 0.061572814777)),
               ('H', (2.06801554556, -0.080044835737, 0.0393130891766)))
    R_MGEOS = (R1_MGEO, R2_MGEO)
    P_MGEOS = (P1_MGEO, P2_MGEO)
    print abstractor_xyz_string(R1_MGEO, 4)
    print abstractee_xyz_string(R1_MGEO, 5)
