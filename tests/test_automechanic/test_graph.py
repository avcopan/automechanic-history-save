""" test the automechanic.graph module
"""
from itertools import permutations
import numpy.random
from automechanic import graph


def test__atoms():
    """ test graph.atoms()
    """
    mgrph = (('C', 'C', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                        (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                        (frozenset([0, 1]), 2)]))
    assert graph.atoms(mgrph) == ('C', 'C', 'H', 'H', 'H', 'H')


def test__bonds():
    """ test graph.bonds()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert (graph.bonds(mgrph1)
            == frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                          (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                          (frozenset([0, 1]), 2)]))
    assert graph.bonds(mgrph2) == frozenset()


def test__bond_keys():
    """ test graph.bond_keys()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert (set(graph.bond_keys(mgrph1))
            == set([frozenset([0, 3]), frozenset([1, 4]), frozenset([0, 2]),
                    frozenset([0, 1]), frozenset([1, 5])]))
    assert graph.bond_keys(mgrph2) == ()


def test__bond_orders():
    """ test graph.bond_orders()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert (graph.bond_orders(mgrph1)
            == {frozenset([1, 5]): 1, frozenset([1, 4]): 1,
                frozenset([0, 2]): 1, frozenset([0, 3]): 1,
                frozenset([0, 1]): 2})
    assert graph.bond_orders(mgrph2) == {}


def test__radical_sites():
    """ test graph.radical_sites()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    mgrph3 = (('C', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([8, 3]), 1),
                         (frozenset([4, 6]), 1), (frozenset([4, 5]), 1),
                         (frozenset([0, 2]), 1), (frozenset([3, 4]), 1),
                         (frozenset([3, 7]), 1), (frozenset([0, 1]), 1)]))
    assert graph.radical_sites(mgrph1) == ()
    assert graph.radical_sites(mgrph2) == (0,)
    assert graph.radical_sites(mgrph3) == (0, 4)


def test__atom_bonds():
    """ test graph.atom_bonds()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert (graph.atom_bonds(mgrph1, 0)
            == frozenset([(frozenset([0, 3]), 1), (frozenset([0, 2]), 1),
                          (frozenset([0, 1]), 2)]))
    assert graph.atom_bonds(mgrph2, 0) == frozenset()


def test__delete_atom():
    """ test graph.delete_atom()
    """
    mgrph = (('C', 'O', 'H', 'H'),
             frozenset([(frozenset([0, 3]), 1), (frozenset([0, 2]), 1),
                        (frozenset([0, 1]), 2)]))
    assert graph.delete_atom(mgrph, 0) == (('O', 'H', 'H'), frozenset([]))
    assert graph.delete_atom(mgrph, 1) == (('C', 'H', 'H'),
                                           frozenset([(frozenset([0, 2]), 1),
                                                      (frozenset([0, 1]), 1)]))
    assert graph.delete_atom(mgrph, 2) == (('C', 'O', 'H'),
                                           frozenset([(frozenset([0, 2]), 1),
                                                      (frozenset([0, 1]), 2)]))


def test__bind_atom():
    """ test graph.bind_atom()
    """
    mgrph = (('O', 'H'), frozenset([(frozenset([0, 1]), 1)]))
    assert (graph.bind_atom(mgrph, 0, 'H', order=1)
            == (('O', 'H', 'H'), frozenset([(frozenset([0, 2]), 1),
                                            (frozenset([0, 1]), 1)])))


def test__atom_neighborhood_indices():
    """ test graph.atom_neighborhood_indices()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert graph.atom_neighborhood_indices(mgrph1, 0) == frozenset([1, 2, 3])
    assert graph.atom_neighborhood_indices(mgrph2, 0) == frozenset()


def test__atom_bond_count():
    """ test graph.atom_bond_count()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    assert graph.atom_bond_count(mgrph1, 0) == 4
    assert graph.atom_bond_count(mgrph2, 0) == 0


def test__atom_free_electrons():
    """ test graph.atom_free_electrons()
    """
    mgrph1 = (('C', 'C', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                         (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 1]), 2)]))
    mgrph2 = ('H', frozenset())
    mgrph3 = (('C', 'H', 'H'),
              frozenset([(frozenset([0, 2]), 1), (frozenset([0, 1]), 1)]))
    assert graph.atom_free_electrons(mgrph1, 0) == 0
    assert graph.atom_free_electrons(mgrph2, 0) == 1
    assert graph.atom_free_electrons(mgrph3, 0) == 2


def test__disjoint_union():
    """ test graph.disjoint_union()
    """
    mgrph1 = (('O',), frozenset())
    mgrph2 = (('C', 'H', 'H'),
              frozenset([(frozenset([0, 2]), 1), (frozenset([0, 1]), 1)]))
    mgrph3 = (('H',), frozenset())
    assert (graph.disjoint_union(mgrph1, mgrph2)
            == (('O', 'C', 'H', 'H'),
                frozenset([(frozenset([1, 3]), 1), (frozenset([1, 2]), 1)])))
    assert (graph.disjoint_union(mgrph2, mgrph3)
            == (('C', 'H', 'H', 'H'),
                frozenset([(frozenset([0, 2]), 1), (frozenset([0, 1]), 1)])))


def test__bind():
    """ test graph.bind()
    """
    mgrph1 = (('O',), frozenset())
    mgrph2 = (('C', 'H', 'H'),
              frozenset([(frozenset([0, 2]), 1), (frozenset([0, 1]), 1)]))
    mgrph3 = (('H',), frozenset())
    assert (graph.bind(mgrph1, mgrph2, 0, 0, order=2)
            == (('O', 'C', 'H', 'H'),
                frozenset([(frozenset([1, 3]), 1), (frozenset([1, 2]), 1),
                           (frozenset([0, 1]), 2)])))
    assert (graph.bind(mgrph2, mgrph3, 0, 0, order=1)
            == (('C', 'H', 'H', 'H'),
                frozenset([(frozenset([0, 3]), 1), (frozenset([0, 2]), 1),
                           (frozenset([0, 1]), 1)])))


def test__permute_atoms():
    """ test graph.permute_atoms()
    """
    mgrph = (('C', 'C', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                        (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                        (frozenset([0, 1]), 2)]))
    assert (graph.permute_atoms(mgrph, (5, 2, 1, 4, 3, 0)) ==
            (('H', 'H', 'C', 'H', 'H', 'C'),
             frozenset([(frozenset([2, 5]), 2), (frozenset([2, 3]), 1),
                        (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                        (frozenset([4, 5]), 1)])))


def test__isomorphic():
    """ test graph.isomorphic()
    """
    mgrph = (('C', 'C', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([0, 3]), 1), (frozenset([1, 4]), 1),
                        (frozenset([0, 2]), 1), (frozenset([1, 5]), 1),
                        (frozenset([0, 1]), 2)]))
    assert all(graph.isomorphic(mgrph, graph.permute_atoms(mgrph, p))
               for p in permutations(range(6)))


def test__isomorphism():
    """ test graph.isomorphism()
    """
    mgrph = (('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([2, 4]), 1), (frozenset([0, 6]), 1),
                        (frozenset([0, 5]), 1), (frozenset([1, 2]), 1),
                        (frozenset([0, 1]), 2), (frozenset([2, 3]), 1),
                        (frozenset([1, 7]), 1)]))

    natms = len(graph.atoms(mgrph))
    for _ in range(10):
        pmt = tuple(numpy.random.permutation(natms))
        pmt_mgrph = graph.permute_atoms(mgrph, pmt)
        iso = graph.isomorphism(mgrph, pmt_mgrph)
        assert graph.permute_atoms(mgrph, iso) == pmt_mgrph


if __name__ == '__main__':
    test__isomorphism()
    test__isomorphic()
