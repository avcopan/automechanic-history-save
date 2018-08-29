""" a library of specialized geom functions
"""
from itertools import product
from itertools import chain
from .geom import resonance_graphs
from .graph import multibond_opening_resonances
from .graphlib import addition_indices as graph_addition_indices
from .graphlib import (forward_abstraction_indices as
                       graph_forward_abstraction_indices)


def addition_indices(x_mgeo, y_mgeo, xy_mgeo):
    """ find addition for a pre-sorted reaction
    """
    idxs = None

    x_mgrphs = resonance_graphs(x_mgeo)
    y_mgrphs = resonance_graphs(y_mgeo)
    xy_mgrphs = resonance_graphs(xy_mgeo)

    for x_mgrph, y_mgrph, xy_mgrph in product(x_mgrphs, y_mgrphs, xy_mgrphs):
        idxs = graph_addition_indices(x_mgrph, y_mgrph, xy_mgrph)
        if idxs:
            break

    return idxs


def abstraction_indices(r1h_mgeo, r2_mgeo, p1_mgeo, p2h_mgeo):
    """ find hydrogen abstraction for a pre-sorted reaction
    """
    idxs = None
    idxs1 = forward_abstraction_indices(r1h_mgeo, p1_mgeo)
    idxs2 = forward_abstraction_indices(p2h_mgeo, r2_mgeo)
    if idxs1 and idxs2:
        q1h_idx, q1_idx = idxs1
        q2h_idx, q2_idx = idxs2
        idxs = (q1h_idx, q2_idx, q1_idx, q2h_idx)
    return idxs


def forward_abstraction_indices(qh_mgeo, q_mgeo):
    """ identify abstraction indices
    """
    idxs = None

    qh_mgrphs_iter = chain(*map(multibond_opening_resonances,
                                resonance_graphs(qh_mgeo)))
    q_mgrphs_iter = chain(*map(multibond_opening_resonances,
                               resonance_graphs(q_mgeo)))

    for qh_mgrph, q_mgrph in product(qh_mgrphs_iter, q_mgrphs_iter):
        idxs = graph_forward_abstraction_indices(qh_mgrph, q_mgrph)
        if idxs:
            break

    return idxs
