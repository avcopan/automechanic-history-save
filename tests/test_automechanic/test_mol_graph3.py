""" test the automechanc.mol.graph module
"""
from automechanic.mol import graph3 as graph


C8H13O_CGR = (
    {0: ('C', 3), 1: ('C', 3), 2: ('C', 1), 3: ('C', 1), 4: ('C', 1),
     5: ('C', 1), 6: ('C', 2), 7: ('C', 1), 8: ('O', 0)},
    {frozenset({0, 2}): (1,), frozenset({1, 3}): (1,),
     frozenset({2, 4}): (1,), frozenset({3, 5}): (1,),
     frozenset({4, 6}): (1,), frozenset({5, 7}): (1,),
     frozenset({6, 7}): (1,), frozenset({8, 7}): (1,)})
C8H13O_RGR = (
    {0: ('C', 3), 1: ('C', 3), 2: ('C', 1), 3: ('C', 1), 4: ('C', 1),
     5: ('C', 1), 6: ('C', 2), 7: ('C', 1), 8: ('O', 0)},
    {frozenset({0, 2}): (1,), frozenset({1, 3}): (1,),
     frozenset({2, 4}): (2,), frozenset({3, 5}): (2,),
     frozenset({4, 6}): (1,), frozenset({5, 7}): (1,),
     frozenset({6, 7}): (1,), frozenset({8, 7}): (1,)})
C8H13O_SGR = (
    {0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 2, None), 7: ('C', 1, False), 8: ('O', 0, None)},
    {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
     frozenset({2, 4}): (1, False), frozenset({3, 5}): (1, False),
     frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None)})
CH2FH2H_CGR_IMP = (
    {1: ('F', 0), 3: ('C', 2), 4: ('H', 1), 6: ('H', 0)},
    {frozenset({1, 3}): (1,)})
CH2FH2H_CGR_EXP = (
    {0: ('H', 0), 1: ('F', 0), 2: ('H', 0), 3: ('C', 0),
     4: ('H', 0), 5: ('H', 0), 6: ('H', 0)},
    {frozenset({1, 3}): (1,), frozenset({2, 3}): (1,),
     frozenset({0, 3}): (1,), frozenset({4, 5}): (1,)})


def test__from_data():
    """ test graph.from_data

    also tests a bunch of accessors
    """
    assert graph.from_data(
        graph.atom_keys(CH2FH2H_CGR_EXP),
        graph.bond_keys(CH2FH2H_CGR_EXP),
        graph.atom_symbols(CH2FH2H_CGR_EXP)
    ) == CH2FH2H_CGR_EXP
    assert graph.from_data(
        graph.atom_keys(C8H13O_CGR),
        graph.bond_keys(C8H13O_CGR),
        graph.atom_symbols(C8H13O_CGR),
        atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(C8H13O_CGR)
    ) == C8H13O_CGR
    assert graph.from_data(
        graph.atom_keys(C8H13O_RGR),
        graph.bond_keys(C8H13O_RGR),
        graph.atom_symbols(C8H13O_RGR),
        atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(C8H13O_CGR),
        bnd_ord_dct=graph.res.bond_orders(C8H13O_RGR)
    ) == C8H13O_RGR
    assert graph.from_data(
        graph.atom_keys(C8H13O_SGR),
        graph.bond_keys(C8H13O_SGR),
        graph.atom_symbols(C8H13O_SGR),
        atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(C8H13O_SGR),
        atm_par_dct=graph.stereo.atom_parities(C8H13O_SGR),
        bnd_par_dct=graph.stereo.bond_parities(C8H13O_SGR)
    ) == C8H13O_SGR


def test__set_atom_implicit_hydrogen_valences():
    """ test graph.set_atom_implicit_hydrogen_valences
    """
    assert graph.set_atom_implicit_hydrogen_valences(
        CH2FH2H_CGR_IMP, {3: 1, 6: 1}
    ) == ({1: ('F', 0), 3: ('C', 1), 4: ('H', 1), 6: ('H', 1)},
          {frozenset({1, 3}): (1,)})


def test__atom_nuclear_charges():
    """ test graph.atom_nuclear_charges
    """
    assert (graph.atom_nuclear_charges(C8H13O_CGR) ==
            {0: 6, 1: 6, 2: 6, 3: 6, 4: 6, 5: 6, 6: 6, 7: 6, 8: 8})


def test__atom_total_valences():
    """ test graph.atom_total_valences
    """
    assert (graph.atom_total_valences(C8H13O_CGR) ==
            {0: 4, 1: 4, 2: 4, 3: 4, 4: 4, 5: 4, 6: 4, 7: 4, 8: 2})


def test__atom_bonds():
    """ test graph.atom_bonds
    """
    assert (graph.atom_bonds(C8H13O_CGR) ==
            {0: {frozenset({0, 2}): (1,)},
             1: {frozenset({1, 3}): (1,)},
             2: {frozenset({0, 2}): (1,), frozenset({2, 4}): (1,)},
             3: {frozenset({1, 3}): (1,), frozenset({3, 5}): (1,)},
             4: {frozenset({2, 4}): (1,), frozenset({4, 6}): (1,)},
             5: {frozenset({3, 5}): (1,), frozenset({5, 7}): (1,)},
             6: {frozenset({4, 6}): (1,), frozenset({6, 7}): (1,)},
             7: {frozenset({5, 7}): (1,), frozenset({6, 7}): (1,),
                 frozenset({8, 7}): (1,)},
             8: {frozenset({8, 7}): (1,)}})


def test__atom_neighbor_keys():
    """ test graph.atom_neighbor_keys
    """
    assert (graph.atom_neighbor_keys(C8H13O_CGR) ==
            {0: (2,), 1: (3,), 2: (0, 4), 3: (1, 5), 4: (2, 6), 5: (3, 7),
             6: (4, 7), 7: (5, 6, 8), 8: (7,)})


# change this to `atom_explicit_hydrogen_valences`
# with separate `backbone_keys` and `explicit_hydrogen_keys` functions
# def test__atom_explicit_hydrogen_keys():
#     """ test graph.atom_explicit_hydrogen_keys
#     """
#     assert (graph.atom_explicit_hydrogen_keys(CH2FH2H_CGR_EXP) ==
#             {3: (0, 2), 4: (5,), 5: (4,)})


def test__res__no_pi():
    """ test graph.res.no_pi
    """
    no_pi_rgr = graph.res.set_bond_orders(
        C8H13O_RGR, {frozenset({2, 4}): 1, frozenset({3, 5}): 1})
    assert graph.res.no_pi(C8H13O_RGR) == C8H13O_CGR == no_pi_rgr


def test__stereo__no_centers():
    """ test graph.stereo.no_centers
    """
    no_ctr_sgr = graph.stereo.set_bond_parities(
        graph.stereo.set_atom_parities(C8H13O_SGR, {7: None}),
        {frozenset({2, 4}): None, frozenset({3, 5}): None})
    assert graph.stereo.no_centers(C8H13O_SGR) == no_ctr_sgr


if __name__ == '__main__':
    test__from_data()
    test__set_atom_implicit_hydrogen_valences()
    test__atom_nuclear_charges()
    test__atom_total_valences()
    test__atom_bonds()
    test__atom_neighbor_keys()
    test__res__no_pi()
    test__stereo__no_centers()
