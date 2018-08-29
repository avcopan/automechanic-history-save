""" a library of specialized geom functions
"""
from itertools import product
from itertools import permutations
from .graph import atoms
from .graph import bind
from .graph import bind_atom
from .graph import isomorphism
from .graph import radical_sites
from .graph import multibond_keys
from .graph import increment_bond_order


def forward_abstraction_indices(qh_mgrph, q_mgrph):
    """ find abstraction indices
    """
    idxs = None
    for idx in radical_sites(q_mgrph):
        qh_mgrph_ = bind_atom(q_mgrph, idx, 'H')
        iso = isomorphism(qh_mgrph, qh_mgrph_)
        if iso:
            qh_idx = int(iso[-1])
            q_idx = int(idx)
            idxs = (qh_idx, q_idx)
    return idxs


def addition_indices(x_mgrph, y_mgrph, xy_mgrph):
    """ find addition indices
    """
    idxs = None

    x_idxs = radical_sites(x_mgrph)
    y_bkeys = multibond_keys(y_mgrph)
    for x_idx, y_bkey in product(x_idxs, y_bkeys):
        for y_idx, y_idx_other in permutations(y_bkey):
            y_open_mgrph = increment_bond_order(
                y_mgrph, y_idx, y_idx_other, incr=-1)
            xy_mgrph_ = bind(x_mgrph, y_open_mgrph, x_idx, y_idx)
            iso = isomorphism(xy_mgrph, xy_mgrph_)
            if iso:
                natms_x = len(atoms(x_mgrph))
                xy_idx1 = int(iso[x_idx])
                xy_idx2 = int(iso[natms_x + y_idx])
                idxs = (x_idx, y_idx, xy_idx1, xy_idx2)

    return idxs
