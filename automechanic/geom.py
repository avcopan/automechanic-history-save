""" cartesian geometry functions
"""
from itertools import product
from itertools import chain
from .ipybel.geom import smiles
from .ipyx2z.geom import graph
from .ipyx2z.geom import resonance_graphs
from .ipyx2z.geom import radical_sites
from .graph import multibond_opening_resonances
from .graph import addition_indices as graph_addition_indices
from .graph import migration_indices as graph_migration_indices
from .graph import (forward_abstraction_indices as
                    graph_forward_abstraction_indices)


def xyz_string(mgeo, labels=None):
    """ .xyz format string for this cartesian geometry

    :param labels: optional labels for the beginnings of atom lines, by index
    :type labels: dict
    """
    labels = {} if labels is None else labels

    natms = len(mgeo)
    dxyz = '{:d}\n\n'.format(natms)
    for idx, (asymb, xyz) in enumerate(mgeo):
        if idx in labels:
            dxyz += '{:s} '.format(repr(labels[idx]))
        dxyz += '{:s} {:s} {:s} {:s}\n'.format(asymb, *map(repr, xyz))
    return dxyz


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


def migration_indices(r_mgeo, p_mgeo):
    """ identify migration indices
    """
    idxs = None

    r_mgrphs_iter = chain(*map(multibond_opening_resonances,
                               resonance_graphs(r_mgeo)))
    p_mgrphs_iter = chain(*map(multibond_opening_resonances,
                               resonance_graphs(p_mgeo)))

    for r_mgrph, p_mgrph in product(r_mgrphs_iter, p_mgrphs_iter):
        idxs = graph_migration_indices(r_mgrph, p_mgrph)
        if idxs:
            break

    return idxs


__all__ = ['smiles', 'xyz_string', 'graph', 'resonance_graphs',
           'radical_sites']
