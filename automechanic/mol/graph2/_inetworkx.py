""" networkx interface
"""
import networkx
from ._perm import inverse as _inverse_permutation


def from_automechanic_graph(gra):
    """ newtorkx graph object from an automechanic graph object
    """
    vtcs, edgs = gra
    vtcs = dict(enumerate(vtcs))

    nxg = networkx.Graph()
    for vkey, vprops in vtcs.items():
        nxg.add_node(vkey, props=vprops)
    for ekey, eprops in edgs.items():
        nxg.add_edge(*ekey, props=eprops)
    return nxg


def isomorphism(nxg1, nxg2):
    """ graph isomorphism
    """

    def _same_props(dct1, dct2):
        return dct1['props'] == dct2['props']

    matcher = networkx.algorithms.isomorphism.GraphMatcher(
        nxg1, nxg2, node_match=_same_props, edge_match=_same_props)

    iso = None
    if matcher.is_isomorphic():
        iso_dct = matcher.mapping
        iso_inv = tuple(map(iso_dct.__getitem__, range(len(iso_dct))))
        iso = _inverse_permutation(iso_inv)

    return iso
