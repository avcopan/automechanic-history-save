""" test the automechanic.mol.graph module
"""
from automechanic import mol
from automechanic.mol import graph2 as graph
from automechanic.mol import graph as old_graph


def test__delete_vertices():
    """ test graph.delete_vertices
    """
    gra = ((('H', 0), ('F', 0), ('H', 0), ('C', 0), ('H', 0), ('H', 0)),
           {frozenset({1, 3}): None, frozenset({2, 3}): None,
            frozenset({0, 3}): None, frozenset({4, 5}): None})
    assert (graph.delete_vertices(gra, (3,))
            == ((('H', 0), ('F', 0), ('H', 0), ('H', 0), ('H', 0)),
                {frozenset({3, 4}): None}))


def test__vertex_neighbor_keys():
    """ test graph.vertex_neighbor_keys
    """
    gra = ((('H', 0), ('F', 0), ('H', 0), ('C', 0)),
           {frozenset({1, 3}): None, frozenset({2, 3}): None,
            frozenset({0, 3}): None})
    assert graph.vertex_neighbor_keys(gra, 3) == (0, 1, 2)


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
    ch2f1_cgr = ((('C', 2), ('F', 0)), {frozenset([0, 1]): None})
    ch2f1_mlf = graph.conn.molfile(ch2f1_cgr)
    ch2f1_ich = mol.molfile.inchi(ch2f1_mlf)
    ch2f1_cgr2 = mol.inchi.connectivity_graph(ch2f1_ich)
    print(ch2f1_mlf)
    print(ch2f1_ich)
    print(ch2f1_cgr2)


def test__old_res__radical_sites():
    """ test res.radical_sites()
    """
    rgr = (('C', 'C', 'H', 'H'),
           {frozenset([0, 3]): 1, frozenset([0, 2]): 1, frozenset([0, 1]): 2})
    assert old_graph.res.radical_sites(rgr) == {1: 2}


def test__old_conn__possible_spin_multiplicities():
    """ test conn.possible_spin_multiplicities
    """
    ch0_cgr = (('C'), {})
    ch1_cgr = (('C', 'H'), {frozenset([0, 1]): None})
    ch2_cgr = (('C', 'H', 'H'),
               {frozenset([0, 1]): None, frozenset([0, 2]): None})
    ch3_cgr = (('C', 'H', 'H', 'H'),
               {frozenset([0, 1]): None, frozenset([0, 2]): None,
                frozenset([0, 3]): None})
    ch4_cgr = (('C', 'H', 'H', 'H', 'H'),
               {frozenset([0, 1]): None, frozenset([0, 2]): None,
                frozenset([0, 3]): None, frozenset([0, 4]): None})
    o2_cgr = (('O', 'O'), {frozenset([0, 1]): None})
    assert old_graph.conn.possible_spin_multiplicities(ch0_cgr) == (1, 3, 5)
    assert old_graph.conn.possible_spin_multiplicities(ch1_cgr) == (2, 4)
    assert old_graph.conn.possible_spin_multiplicities(ch2_cgr) == (1, 3)
    assert old_graph.conn.possible_spin_multiplicities(ch3_cgr) == (2,)
    assert old_graph.conn.possible_spin_multiplicities(ch4_cgr) == (1,)
    assert old_graph.conn.possible_spin_multiplicities(o2_cgr) == (1, 3)


if __name__ == '__main__':
    test__old_res__radical_sites()
    test__old_conn__possible_spin_multiplicities()
    test__conn__possible_spin_multiplicities()
    # test__conn__molfile()
    test__conn__make_hydrogens_implicit()
    test__permute_vertices()
    test__vertex_neighbor_keys()
    test__delete_vertices()
