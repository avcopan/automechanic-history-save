""" test the automechanic.mol.graph module
"""
from automechanic.mol import graph


def test__graph__radical_sites():
    """ test graph.res.radical_sites()
    """
    rgr = (('C', 'C', 'H', 'H'),
           {frozenset([0, 3]): 1, frozenset([0, 2]): 1, frozenset([0, 1]): 2})
    assert graph.res.radical_sites(rgr) == {1: 2}


def test__connn__possible_spin_multiplicities():
    """ test graph.possible_spin_multiplicities
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
    assert graph.conn.possible_spin_multiplicities(ch0_cgr) == (1, 3, 5)
    assert graph.conn.possible_spin_multiplicities(ch1_cgr) == (2, 4)
    assert graph.conn.possible_spin_multiplicities(ch2_cgr) == (1, 3)
    assert graph.conn.possible_spin_multiplicities(ch3_cgr) == (2,)
    assert graph.conn.possible_spin_multiplicities(ch4_cgr) == (1,)
    assert graph.conn.possible_spin_multiplicities(o2_cgr) == (1, 3)


if __name__ == '__main__':
    test__graph__radical_sites()
    test__connn__possible_spin_multiplicities()
