""" a library of specialized geom functions
"""
from itertools import product
from itertools import chain
from .geom import resonance_graphs
from .graph import multibond_opening_resonances
from .graphlib import abstraction_indices as graph_abstraction_indices


def find_abstraction(r1h_mgeo, r2_mgeo, p1_mgeo, p2h_mgeo):
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

    qh_mgrphs_iter = chain(*map(multibond_opening_resonances,
                                resonance_graphs(qh_mgeo)))
    q_mgrphs_iter = chain(*map(multibond_opening_resonances,
                               resonance_graphs(q_mgeo)))

    for qh_mgrph, q_mgrph in product(qh_mgrphs_iter, q_mgrphs_iter):
        idxs = graph_abstraction_indices(qh_mgrph, q_mgrph)
        if idxs:
            break
    return idxs


if __name__ == '__main__':
    Q1H_MGEO = (('O', (1.160649880891, 0.0627978268174, -0.04083307250753)),
                ('O', (0.670030192617, -0.810012591664, 1.01847301649)),
                ('C', (-0.762764327328, -0.771623494653, 0.971879813498)),
                ('C', (-1.286027334011, 0.630975793141, 1.25724343361)),
                ('C', (-1.286025794377, -1.319989887595, -0.350244634965)),
                ('H', (2.10820223046, -0.0822581224768, 0.135216457591)),
                ('H', (-1.10626335348, -1.433832526376, 1.775585055636)),
                ('H', (-0.959411887323, 1.34277748802, 0.4919543629624)),
                ('H', (-2.37990142664, 0.642772565761, 1.288621209488)),
                ('H', (-0.902576311752, 0.995523287157, 2.215796882617)),
                ('H', (-0.9025737237, -2.330572084284, -0.524722174920)),
                ('H', (-2.379899853887, -1.353045699721, -0.3558229455)),
                ('H', (-0.959410238052, -0.70497322095, -1.195279213213)))
    Q2_MGEO = (('O', (1.07394324121, 0.0501743774966, 0.0668114667428)),
               ('O', (2.355943241199, 0.0501743774966, 0.0668114667428)))
    Q1_MGEO = (('O', (1.01849391329, 0.02212499683403, 0.0463288404118)),
               ('O', (2.424344519778, -0.02924375645281, 0.05658010669234)),
               ('C', (2.80314339128, -1.38272890920, 0.3266805371647)),
               ('C', (2.338504665743, -1.833731474576, 1.707956931019)),
               ('C', (2.338505036357, -2.329364870800, -0.775684415566)),
               ('H', (3.89891588452, -1.38511658503, 0.327157185309)),
               ('H', (1.246260227299, -1.860092545692, 1.776726361401)),
               ('H', (2.71835895531, -2.83374919094, 1.940450620285)),
               ('H', (2.690750920592, -1.13759031891, 2.476303244045)),
               ('H', (2.690751556464, -1.981462639680, -1.752379053941)),
               ('H', (2.7183593275, -3.34202255559, -0.606530163155)),
               ('H', (1.24626061603, -2.380102757813, -0.829068367911)))
    Q2H_MGEO = (('O', (1.09332243633, -0.0516790137299, -0.01278085883290)),
                ('O', (0.596198306921, -1.092054935846, -0.810333075830)),
                ('H', (2.065302405982, -0.05663139196237, -0.01657735215297)))
    print find_abstraction(Q1H_MGEO, Q2_MGEO, Q1_MGEO, Q2H_MGEO)
