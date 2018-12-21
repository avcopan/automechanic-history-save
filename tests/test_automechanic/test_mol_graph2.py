""" test the automechanc.mol.graph module
"""
from automechanic.mol import graph2 as graph


def test__conn__from_data():
    """ test graph.conn.from_data
    """
    cgr = ({0: ('C', 3), 1: ('C', 3), 2: ('C', 1), 3: ('C', 1), 4: ('C', 1),
            5: ('C', 1), 6: ('C', 2), 7: ('C', 1), 8: ('O', 0)},
           {frozenset({0, 2}): None, frozenset({1, 3}): None,
            frozenset({2, 4}): None, frozenset({3, 5}): None,
            frozenset({4, 6}): None, frozenset({5, 7}): None,
            frozenset({6, 7}): None, frozenset({8, 7}): None})
    assert graph.conn.from_data(
        atm_keys=graph.conn.atom_keys(cgr),
        atm_syms=graph.conn.atom_symbols(cgr),
        atm_hyd_cnts=graph.conn.atom_hydrogen_counts(cgr),
        bnd_keys=graph.conn.bond_keys(cgr)
    ) == cgr


def test__res__from_data():
    """ test graph.res.from_data
    """
    rgr = ({0: ('C', 3), 1: ('C', 3), 2: ('C', 1), 3: ('C', 1), 4: ('C', 1),
            5: ('C', 1), 6: ('C', 2), 7: ('C', 1), 8: ('O', 0)},
           {frozenset({0, 2}): 1, frozenset({1, 3}): 1,
            frozenset({2, 4}): 2, frozenset({3, 5}): 2,
            frozenset({4, 6}): 1, frozenset({5, 7}): 1,
            frozenset({6, 7}): 1, frozenset({8, 7}): 1})
    assert graph.res.from_data(
        atm_keys=graph.res.atom_keys(rgr),
        atm_syms=graph.res.atom_symbols(rgr),
        atm_hyd_cnts=graph.res.atom_hydrogen_counts(rgr),
        bnd_keys=graph.res.bond_keys(rgr),
        bnd_ords=graph.res.bond_orders(rgr)
    ) == rgr


def test__stereo__from_data():
    """ test graph.conn.from_data
    """
    sgr = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
            3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
            6: ('C', 2, None), 7: ('C', 1, False), 8: ('O', 0, None)},
           {frozenset({0, 2}): None, frozenset({1, 3}): None,
            frozenset({2, 4}): False, frozenset({3, 5}): False,
            frozenset({4, 6}): None, frozenset({5, 7}): None,
            frozenset({6, 7}): None, frozenset({8, 7}): None})
    assert graph.stereo.from_data(
        atm_keys=graph.stereo.atom_keys(sgr),
        atm_syms=graph.stereo.atom_symbols(sgr),
        atm_hyd_cnts=graph.stereo.atom_hydrogen_counts(sgr),
        atm_pars=graph.stereo.atom_stereo_parities(sgr),
        bnd_keys=graph.stereo.bond_keys(sgr),
        bnd_pars=graph.stereo.bond_stereo_parities(sgr)
    ) == sgr


if __name__ == '__main__':
    test__conn__from_data()
    test__res__from_data()
    test__stereo__from_data()
