""" a library of specialized geom functions
"""
from .graph import bind_atom
from .graph import isomorphism
from .graph import radical_sites


def abstraction_indices(qh_mgrph, q_mgrph):
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
