""" test the automechanc.mol.graph module
"""
import itertools
import numpy
from automechanic.mol import graph3 as graph
import automechanic.mol.graph3.int_xyz_rot as int_xyz_rot


C8H13O_CGR = (
    {0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 2, None), 7: ('C', 1, None), 8: ('O', 0, None)},
    {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
     frozenset({2, 4}): (1, None), frozenset({3, 5}): (1, None),
     frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None)})
C8H13O_RGR = (
    {0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 2, None), 7: ('C', 1, None), 8: ('O', 0, None)},
    {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
     frozenset({2, 4}): (2, None), frozenset({3, 5}): (2, None),
     frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None)})
C8H13O_SGR = (
    {0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
     3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
     6: ('C', 2, None), 7: ('C', 1, False), 8: ('O', 0, None)},
    {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
     frozenset({2, 4}): (1, False), frozenset({3, 5}): (1, False),
     frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
     frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None)})
CH2FH2H_CGR_IMP = (
    {1: ('F', 0, None), 3: ('C', 2, None), 4: ('H', 1, None),
     6: ('H', 0, None)},
    {frozenset({1, 3}): (1, None)})
CH2FH2H_CGR_EXP = (
    {0: ('H', 0, None), 1: ('F', 0, None), 2: ('H', 0, None),
     3: ('C', 0, None), 4: ('H', 0, None), 5: ('H', 0, None),
     6: ('H', 0, None)},
    {frozenset({1, 3}): (1, None), frozenset({2, 3}): (1, None),
     frozenset({0, 3}): (1, None), frozenset({4, 5}): (1, None)})


# test constructors and value getters
def test__from_data():
    """ test graph.from_data

    also tests a bunch of accessors
    """
    assert graph.from_data(
        graph.atom_symbols(CH2FH2H_CGR_EXP),
        graph.bond_keys(CH2FH2H_CGR_EXP)
    ) == CH2FH2H_CGR_EXP
    assert graph.from_data(
        graph.atom_symbols(C8H13O_CGR),
        graph.bond_keys(C8H13O_CGR),
        atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(C8H13O_CGR)
    ) == C8H13O_CGR
    assert graph.from_data(
        graph.atom_symbols(C8H13O_RGR),
        graph.bond_keys(C8H13O_RGR),
        atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(C8H13O_CGR),
        bnd_ord_dct=graph.bond_orders(C8H13O_RGR)
    ) == C8H13O_RGR
    assert graph.from_data(
        graph.atom_symbols(C8H13O_SGR),
        graph.bond_keys(C8H13O_SGR),
        atm_imp_hyd_vlc_dct=graph.atom_implicit_hydrogen_valences(C8H13O_SGR),
        atm_ste_par_dct=graph.atom_stereo_parities(C8H13O_SGR),
        bnd_ste_par_dct=graph.bond_stereo_parities(C8H13O_SGR)
    ) == C8H13O_SGR


# test value setters
def test__set_atom_implicit_hydrogen_valences():
    """ test graph.set_atom_implicit_hydrogen_valences
    """
    assert graph.set_atom_implicit_hydrogen_valences(
        CH2FH2H_CGR_IMP, {3: 1, 4: 0, 6: 1}
    ) == ({1: ('F', 0, None), 3: ('C', 1, None), 4: ('H', 0, None),
           6: ('H', 1, None)},
          {frozenset({1, 3}): (1, None)})


def test__set_atom_stereo_parities():
    """ test graph.set_atom_stereo_parities
    """
    assert graph.atom_stereo_parities(
        graph.set_atom_stereo_parities(C8H13O_CGR, {7: False})
    ) == graph.atom_stereo_parities(C8H13O_SGR)


def test__set_bond_orders():
    """ test graph.set_bond_orders
    """
    assert graph.set_bond_orders(
        C8H13O_CGR, {frozenset({2, 4}): 2, frozenset({3, 5}): 2},
    ) == C8H13O_RGR


def test__set_bond_stereo_parities():
    """ test graph.set_bond_stereo_parities
    """
    assert graph.bond_stereo_parities(
        graph.set_bond_stereo_parities(
            C8H13O_CGR, {frozenset({2, 4}): False, frozenset({3, 5}): False},
        )
    ) == graph.bond_stereo_parities(C8H13O_SGR)


# test derived values
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
            {0: {frozenset({0, 2}): (1, None)},
             1: {frozenset({1, 3}): (1, None)},
             2: {frozenset({0, 2}): (1, None), frozenset({2, 4}): (1, None)},
             3: {frozenset({1, 3}): (1, None), frozenset({3, 5}): (1, None)},
             4: {frozenset({2, 4}): (1, None), frozenset({4, 6}): (1, None)},
             5: {frozenset({3, 5}): (1, None), frozenset({5, 7}): (1, None)},
             6: {frozenset({4, 6}): (1, None), frozenset({6, 7}): (1, None)},
             7: {frozenset({5, 7}): (1, None), frozenset({6, 7}): (1, None),
                 frozenset({8, 7}): (1, None)},
             8: {frozenset({8, 7}): (1, None)}})


def test__atom_neighbor_keys():
    """ test graph.atom_neighbor_keys
    """
    assert (graph.atom_neighbor_keys(C8H13O_CGR) ==
            {0: (2,), 1: (3,), 2: (0, 4), 3: (1, 5), 4: (2, 6), 5: (3, 7),
             6: (4, 7), 7: (5, 6, 8), 8: (7,)})


def test__explicit_hydrogen_keys():
    """ test graph.explicit_hydrogen_keys
    """
    assert graph.explicit_hydrogen_keys(CH2FH2H_CGR_EXP) == (0, 2, 5)


def test__backbone_keys():
    """ test graph.backbone_keys
    """
    assert graph.backbone_keys(CH2FH2H_CGR_EXP) == (1, 3, 4, 6)


def test__atom_explicit_hydrogen_keys():
    """ test graph.atom_explicit_hydrogen_keys
    """
    assert (graph.atom_explicit_hydrogen_keys(CH2FH2H_CGR_EXP) ==
            {0: (), 1: (), 2: (), 3: (0, 2), 4: (5,), 5: (), 6: ()})


def test__ring_keys_list():
    """ test graph.ring_keys_list
    """
    cgr = ({0: ('C', 1, None), 1: ('C', 0, None), 2: ('C', 0, None),
            3: ('C', 0, None), 4: ('C', 0, None), 5: ('N', 2, None),
            6: ('N', 0, None), 7: ('N', 0, None), 8: ('N', 0, None),
            9: ('N', 1, None), 10: ('O', 1, None)},
           {frozenset({10, 4}): (1, None), frozenset({8, 2}): (1, None),
            frozenset({0, 6}): (1, None), frozenset({9, 3}): (1, None),
            frozenset({1, 2}): (1, None), frozenset({3, 7}): (1, None),
            frozenset({2, 5}): (1, None), frozenset({1, 6}): (1, None),
            frozenset({0, 7}): (1, None), frozenset({9, 4}): (1, None),
            frozenset({1, 3}): (1, None), frozenset({8, 4}): (1, None)})
    assert graph.ring_keys_list(cgr) == ((0, 1, 3, 6, 7), (1, 2, 3, 4, 8, 9))


# test transformations
def test__implicit():
    """ test graph.implicit
    """
    assert graph.implicit(CH2FH2H_CGR_EXP) == CH2FH2H_CGR_IMP
    assert graph.implicit(CH2FH2H_CGR_EXP, (1, 3, 4, 6)) == CH2FH2H_CGR_IMP


def test__explicit():
    """ test graph.explicit
    """
    ch2fh2h_cgr_exp = graph.explicit(CH2FH2H_CGR_IMP)
    assert graph.backbone_isomorphic(ch2fh2h_cgr_exp, CH2FH2H_CGR_EXP)
    assert (graph.atom_explicit_hydrogen_keys(ch2fh2h_cgr_exp) ==
            {1: (), 3: (7, 8), 4: (9,), 6: (), 7: (), 8: (), 9: ()})


def test__delete_atoms():
    """ test graph.delete_atoms
    """
    assert (graph.delete_atoms(CH2FH2H_CGR_EXP, (0, 2, 5)) ==
            ({1: ('F', 0, None), 3: ('C', 0, None), 4: ('H', 0, None),
              6: ('H', 0, None)},
             {frozenset({1, 3}): (1, None)}))


def test__add_explicit_hydrogens():
    """ test graph.add_explicit_hydrogens
    """
    assert graph.add_explicit_hydrogens(
        CH2FH2H_CGR_IMP, {3: 2, 4: 1}
    ) == ({1: ('F', 0, None), 3: ('C', 2, None), 4: ('H', 1, None),
           6: ('H', 0, None), 7: ('H', 0, None), 8: ('H', 0, None),
           9: ('H', 0, None)},
          {frozenset({1, 3}): (1, None), frozenset({3, 7}): (1, None),
           frozenset({8, 3}): (1, None), frozenset({9, 4}): (1, None)})


def test__subgraph():
    """ test graph.subgraph
    """
    assert (graph.subgraph(CH2FH2H_CGR_EXP, (1, 3, 4, 6)) ==
            ({1: ('F', 0, None), 3: ('C', 0, None), 4: ('H', 0, None),
              6: ('H', 0, None)},
             {frozenset({1, 3}): (1, None)}))


def test__relabel():
    """ test graph.relabel
    """
    assert graph.relabel(
        CH2FH2H_CGR_IMP, {1: 0, 3: 1, 4: 2, 6: 3}
    ) == ({0: ('F', 0, None), 1: ('C', 2, None), 2: ('H', 1, None),
           3: ('H', 0, None)},
          {frozenset({0, 1}): (1, None)})


# test comparisons
def test__backbone_isomorphic():
    """ test graph.backbone_isomorphic
    """
    assert graph.backbone_isomorphic(CH2FH2H_CGR_EXP, CH2FH2H_CGR_IMP)
    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        assert graph.backbone_isomorphic(cgr, cgr_pmt)


def test__backbone_isomorphism():
    """ test graph.backbone_isomorphism
    """
    assert (graph.backbone_isomorphism(CH2FH2H_CGR_EXP, CH2FH2H_CGR_IMP) ==
            {1: 1, 3: 3, 4: 4, 6: 6})
    cgr = C8H13O_CGR
    natms = len(graph.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.relabel(cgr, pmt_dct)
        assert graph.backbone_isomorphism(cgr, cgr_pmt) == pmt_dct


# test submodules
def test__int_xyz_rot__aligning_rotator():
    """ test int_xyz_rot.aligning_rotator
    """
    for cmp1, cmp2 in itertools.product(range(3), range(3)):
        for val1, val2 in itertools.product([-1, +1], [-1, +1]):
            uint_xyz1 = [0] * 3
            uint_xyz2 = [0] * 3
            uint_xyz1[cmp1] = val1
            uint_xyz2[cmp2] = val2
            uint_xyz1 = tuple(uint_xyz1)
            uint_xyz2 = tuple(uint_xyz2)
            rot_func = int_xyz_rot.aligning_rotator(uint_xyz1, uint_xyz2)
            print(uint_xyz1, uint_xyz2)
            assert rot_func(uint_xyz1) == uint_xyz2


if __name__ == '__main__':
    # test constructors and value getters
    test__from_data()
    # test value setters
    test__set_atom_implicit_hydrogen_valences()
    test__set_atom_stereo_parities()
    test__set_bond_orders()
    test__set_bond_stereo_parities()
    # test derived values
    test__atom_nuclear_charges()
    test__atom_total_valences()
    test__atom_bonds()
    test__atom_neighbor_keys()
    test__explicit_hydrogen_keys()
    test__backbone_keys()
    test__atom_explicit_hydrogen_keys()
    test__ring_keys_list()
    # test transformations
    test__implicit()
    test__explicit()
    test__delete_atoms()
    test__add_explicit_hydrogens()
    test__subgraph()
    test__relabel()
    # test comparisons
    test__backbone_isomorphic()
    test__backbone_isomorphism()
    # test submodules
    test__int_xyz_rot__aligning_rotator()
