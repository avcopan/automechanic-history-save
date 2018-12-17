""" test the automechanic.mol.graph module
"""
import numpy
from automechanic import mol
from automechanic.mol import graph2 as graph


def test__vertex_neighbor_keys():
    """ test graph.vertex_neighbor_keys
    """
    gra = ((('H', 0), ('F', 0), ('H', 0), ('C', 0)),
           {frozenset({1, 3}): None, frozenset({2, 3}): None,
            frozenset({0, 3}): None})
    assert graph.vertex_neighbor_keys(gra, 3) == (0, 1, 2)


def test__isomorphic():
    """ test graph.isomorphic
    """
    catm_cgr = ((('C', 0),), {})
    oatm_cgr = ((('O', 0),), {})
    assert graph.isomorphic(catm_cgr, catm_cgr)
    assert graph.isomorphic(catm_cgr, oatm_cgr) is False


def test__isomorphism():
    """ test graph.isomorphism
    """
    gra = ((('C', 3), ('C', 3), ('C', 1), ('C', 1), ('C', 1), ('C', 1),
            ('C', 2), ('C', 1), ('O', 0)),
           {frozenset({4, 6}): None, frozenset({6, 7}): None,
            frozenset({0, 2}): None, frozenset({8, 7}): None,
            frozenset({2, 4}): None, frozenset({3, 5}): None,
            frozenset({1, 3}): None, frozenset({5, 7}): None})
    nvtcs = len(graph.vertices(gra))
    for _ in range(10):
        pmt = tuple(numpy.random.permutation(nvtcs))
        gra_pmt = graph.permute_vertices(gra, pmt)
        iso = graph.isomorphism(gra, gra_pmt)
        assert graph.permute_vertices(gra, iso) == gra_pmt


def test__induced_subgraph():
    """ test graph.induced_subgraph
    """
    gra = ((('H', 0), ('F', 0), ('H', 0), ('C', 0), ('H', 0), ('H', 0)),
           {frozenset({1, 3}): None, frozenset({2, 3}): None,
            frozenset({0, 3}): None, frozenset({4, 5}): None})
    assert (graph.induced_subgraph(gra, (1, 3, 4, 5))
            == ((('F', 0), ('C', 0), ('H', 0), ('H', 0)),
                {frozenset({2, 3}): None, frozenset({0, 1}): None}))


def test__delete_vertices():
    """ test graph.delete_vertices
    """
    gra = ((('H', 0), ('F', 0), ('H', 0), ('C', 0), ('H', 0), ('H', 0)),
           {frozenset({1, 3}): None, frozenset({2, 3}): None,
            frozenset({0, 3}): None, frozenset({4, 5}): None})
    assert (graph.delete_vertices(gra, (3,))
            == ((('H', 0), ('F', 0), ('H', 0), ('H', 0), ('H', 0)),
                {frozenset({3, 4}): None}))


def test__branch():
    """ test graph.branch
    """
    gra = ((('C', 1), ('C', 1), ('C', 2), ('C', 1), ('F', 0)),
           {frozenset({3, 4}): None, frozenset({2, 3}): None,
            frozenset({0, 1}): None, frozenset({0, 2}): None,
            frozenset({1, 3}): None})
    assert (graph.branch(gra, 3, 1) ==
            ((('C', 1), ('C', 1), ('C', 2), ('C', 1)),
             {frozenset({2, 3}): None, frozenset({0, 1}): None,
              frozenset({0, 2}): None, frozenset({1, 3}): None}))
    assert (graph.branch(gra, 3, 2) ==
            ((('C', 1), ('C', 1), ('C', 2), ('C', 1)),
             {frozenset({2, 3}): None, frozenset({0, 1}): None,
              frozenset({0, 2}): None, frozenset({1, 3}): None}))
    assert (graph.branch(gra, 3, 4) ==
            ((('C', 1), ('F', 0)), {frozenset({0, 1}): None}))


def test__ring_keys_list():
    """ test graph.ring_keys_list
    """
    gra = ((('C', 1), ('C', 0), ('C', 0), ('C', 0), ('C', 0), ('N', 2),
            ('N', 0), ('N', 0), ('N', 0), ('N', 1), ('O', 1)),
           {frozenset({10, 4}): None, frozenset({8, 2}): None,
            frozenset({0, 6}): None, frozenset({9, 3}): None,
            frozenset({1, 2}): None, frozenset({3, 7}): None,
            frozenset({2, 5}): None, frozenset({1, 6}): None,
            frozenset({0, 7}): None, frozenset({9, 4}): None,
            frozenset({1, 3}): None, frozenset({8, 4}): None})
    print(graph.ring_keys_list(gra))


def test__permute_vertices():
    """ test graph.permute_vertices
    """
    gra = ((('H', 0), ('F', 0), ('H', 0), ('C', 0)),
           {frozenset({1, 3}): None, frozenset({2, 3}): None,
            frozenset({0, 3}): None})
    assert (graph.permute_vertices(gra, (3, 1, 2, 0))
            == ((('C', 0), ('F', 0), ('H', 0), ('H', 0)),
                {frozenset({0, 1}): None, frozenset({0, 2}): None,
                 frozenset({0, 3}): None}))


def test__conn__make_hydrogens_implicit():
    """ test graph.conn.make_hydrogens_implicit
    """
    gra = ((('H', 0), ('F', 0), ('H', 0), ('C', 0), ('H', 0), ('H', 0)),
           {frozenset({1, 3}): None, frozenset({2, 3}): None,
            frozenset({0, 3}): None, frozenset({4, 5}): None})
    assert (graph.conn.make_hydrogens_implicit(gra)
            == ((('F', 0), ('C', 2), ('H', 1)), {frozenset({0, 1}): None}))


def test__conn__potential_pi_bond_keys():
    """ test graph.conn.potential_pi_bond_keys
    """
    cgr = ((('C', 3), ('C', 3), ('C', 1), ('C', 1), ('C', 1), ('C', 1),
            ('C', 2), ('C', 1), ('O', 0)),
           {frozenset({4, 6}): None, frozenset({6, 7}): None,
            frozenset({0, 2}): None, frozenset({8, 7}): None,
            frozenset({2, 4}): None, frozenset({3, 5}): None,
            frozenset({1, 3}): None, frozenset({5, 7}): None})
    assert (graph.conn.potential_pi_bond_keys(cgr)
            == (frozenset({2, 4}), frozenset({3, 5})))


def test__conn__stereogenic_atoms():
    """ test graph.conn.stereogenic_atoms
    """
    cgr = ((('C', 3), ('C', 3), ('C', 1), ('C', 1), ('C', 1), ('C', 1),
            ('C', 2), ('C', 1), ('O', 0)),
           {frozenset({4, 6}): None, frozenset({6, 7}): None,
            frozenset({0, 2}): None, frozenset({8, 7}): None,
            frozenset({2, 4}): None, frozenset({3, 5}): None,
            frozenset({1, 3}): None, frozenset({5, 7}): None})
    assert graph.conn.stereogenic_atoms(cgr) == (7,)

    cgr = ((('C', 3), ('C', 1), ('C', 0), ('C', 1), ('F', 0)),
           {frozenset({3, 4}): None, frozenset({2, 3}): None,
            frozenset({1, 2}): None, frozenset({0, 2}): None,
            frozenset({1, 3}): None})
    assert graph.conn.stereogenic_atoms(cgr) == (3,)

    cgr = ((('C', 3), ('C', 1), ('C', 0), ('C', 2)),
           {frozenset({2, 3}): None, frozenset({1, 2}): None,
            frozenset({0, 2}): None, frozenset({1, 3}): None})
    assert graph.conn.stereogenic_atoms(cgr) == ()


def test__conn__stereogenic_bonds():
    """ test graph.conn.stereogenic_bonds
    """
    cgr0 = ((('C', 2), ('C', 0), ('Cl', 0), ('F', 0)),
            {frozenset({0, 1}): None, frozenset({1, 3}): None,
             frozenset({1, 2}): None})
    cgr1 = ((('C', 1), ('C', 0), ('Cl', 0), ('F', 0), ('F', 0)),
            {frozenset({0, 1}): None, frozenset({1, 4}): None,
             frozenset({1, 2}): None, frozenset({0, 3}): None})
    cgr2 = ((('C', 3), ('C', 3), ('C', 1), ('C', 1), ('C', 1), ('C', 1),
             ('C', 2), ('C', 1), ('O', 0)),
            {frozenset({4, 6}): None, frozenset({6, 7}): None,
             frozenset({0, 2}): None, frozenset({8, 7}): None,
             frozenset({2, 4}): None, frozenset({3, 5}): None,
             frozenset({1, 3}): None, frozenset({5, 7}): None})
    assert graph.conn.stereogenic_bonds(cgr0) == ()
    assert graph.conn.stereogenic_bonds(cgr1) == (frozenset({0, 1}),)
    assert graph.conn.stereogenic_bonds(cgr2) == (frozenset({2, 4}),
                                                  frozenset({3, 5}))


def test__conn__possible_spin_multiplicities():
    """ test graph.conn.possible_spin_multiplicities
    """
    catm_cgr = ((('C', 0),), {})
    ch0f1_cgr = ((('C', 0), ('F', 0)), {frozenset([0, 1]): None})
    ch1f1_cgr = ((('C', 1), ('F', 0)), {frozenset([0, 1]): None})
    ch2f1_cgr = ((('C', 2), ('F', 0)), {frozenset([0, 1]): None})
    ch2f2_cgr = ((('C', 2), ('F', 0), ('F', 0)),
                 {frozenset([0, 1]): None, frozenset([0, 2]): None})
    o2_cgr = ((('O', 0), ('O', 0)), {frozenset([0, 1]): None})
    assert graph.conn.possible_spin_multiplicities(catm_cgr) == (1, 3, 5)
    assert graph.conn.possible_spin_multiplicities(ch0f1_cgr) == (2, 4)
    assert graph.conn.possible_spin_multiplicities(ch1f1_cgr) == (1, 3)
    assert graph.conn.possible_spin_multiplicities(ch2f1_cgr) == (2,)
    assert graph.conn.possible_spin_multiplicities(ch2f2_cgr) == (1,)
    assert graph.conn.possible_spin_multiplicities(o2_cgr) == (1, 3)


def test__conn__molfile():
    """ test graph.conn.molfile
    """
    # RDKit can't handle these, so we need a workaround:
    # ch1_cgr = ((('C', 0),), {})
    # ch1_cgr = ((('C', 1),), {})
    ch2_cgr = ((('C', 2),), {})
    ch3_cgr = ((('C', 3),), {})
    ch2f1_cgr = ((('C', 2), ('F', 0)), {frozenset([0, 1]): None})
    c2h2f2_cgr = ((('C', 1), ('C', 1), ('F', 0), ('F', 0)),
                  {frozenset({0, 1}): None, frozenset({0, 2}): None,
                   frozenset({1, 3}): None})
    assert (mol.molfile.inchi(graph.conn.molfile(ch2f1_cgr))
            == 'InChI=1S/CH2F/c1-2/h1H2')
    assert (mol.molfile.inchi(graph.conn.molfile(ch2_cgr))
            == 'InChI=1S/CH2/h1H2')
    assert (mol.molfile.inchi(graph.conn.molfile(ch3_cgr))
            == 'InChI=1S/CH3/h1H3')
    assert (mol.molfile.inchi(graph.conn.molfile(c2h2f2_cgr))
            == 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H')


if __name__ == '__main__':
    test__vertex_neighbor_keys()
    test__isomorphism()
    test__induced_subgraph()
    test__delete_vertices()
    test__branch()
    test__ring_keys_list()
    test__permute_vertices()
    test__conn__possible_spin_multiplicities()
    test__conn__make_hydrogens_implicit()
    test__conn__potential_pi_bond_keys()
    test__conn__stereogenic_atoms()
    test__conn__stereogenic_bonds()
    test__conn__molfile()
