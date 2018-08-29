""" test the automechanic.graphlib module
"""
from automechanic import graph
from automechanic import graphlib


def test__forward_abstraction_indices():
    """ test graphlib.forward_abstraction_indices()
    """
    qh_mgrph = (('C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
                frozenset([(frozenset([0, 6]), 1), (frozenset([1, 2]), 1),
                           (frozenset([0, 5]), 1), (frozenset([8, 1]), 1),
                           (frozenset([1, 7]), 1), (frozenset([2, 3]), 2),
                           (frozenset([10, 3]), 1), (frozenset([0, 4]), 1),
                           (frozenset([3, 11]), 1), (frozenset([9, 2]), 1),
                           (frozenset([0, 1]), 1)]))
    q_mgrph = (('C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
               frozenset([(frozenset([0, 6]), 1), (frozenset([3, 5]), 1),
                          (frozenset([0, 7]), 1), (frozenset([9, 1]), 1),
                          (frozenset([3, 4]), 1), (frozenset([2, 3]), 2),
                          (frozenset([8, 0]), 1), (frozenset([0, 1]), 1),
                          (frozenset([10, 2]), 1), (frozenset([1, 2]), 1)]))
    qh_idx, q_idx = graphlib.forward_abstraction_indices(qh_mgrph, q_mgrph)
    assert (qh_idx, q_idx) == (8, 1)
    assert graph.isomorphism(q_mgrph, graph.delete_atom(qh_mgrph, qh_idx))
    assert graph.isomorphism(qh_mgrph, graph.bind_atom(q_mgrph, q_idx, 'H'))


def test__addition_indices():
    """ test graphlib.addition_indices()
    """
    x_mgrph = (('O', 'H'), frozenset([(frozenset([0, 1]), 1)]))
    y_mgrph = (('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'),
               frozenset([(frozenset([0, 3]), 1), (frozenset([1, 5]), 1),
                          (frozenset([1, 2]), 1), (frozenset([0, 4]), 1),
                          (frozenset([8, 2]), 1), (frozenset([2, 6]), 1),
                          (frozenset([0, 1]), 2), (frozenset([2, 7]), 1)]))
    xy_mgrph = (('C', 'H', 'H', 'C', 'O', 'C', 'H', 'H', 'H', 'H', 'H'),
                frozenset([(frozenset([10, 5]), 1), (frozenset([3, 5]), 1),
                           (frozenset([3, 6]), 1), (frozenset([4, 7]), 1),
                           (frozenset([0, 2]), 1), (frozenset([3, 4]), 1),
                           (frozenset([0, 3]), 1), (frozenset([8, 5]), 1),
                           (frozenset([9, 5]), 1), (frozenset([0, 1]), 1)]))
    idxs = graphlib.addition_indices(x_mgrph, y_mgrph, xy_mgrph)
    assert idxs == (0, 1, 4, 3)


if __name__ == '__main__':
    test__forward_abstraction_indices()
    test__addition_indices()
