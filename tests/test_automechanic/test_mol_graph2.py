""" test the automechanc.mol.graph module
"""
import numpy
from automechanic.mol import graph2 as graph

C8H13O_ICH_NO_STEREO = 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3'
C8H13O_ICH = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-/m0/s1')
C8H13O_CGR = ({0: ('C', 3), 1: ('C', 3), 2: ('C', 1), 3: ('C', 1), 4: ('C', 1),
               5: ('C', 1), 6: ('C', 2), 7: ('C', 1), 8: ('O', 0)},
              {frozenset({0, 2}): None, frozenset({1, 3}): None,
               frozenset({2, 4}): None, frozenset({3, 5}): None,
               frozenset({4, 6}): None, frozenset({5, 7}): None,
               frozenset({6, 7}): None, frozenset({8, 7}): None})
C8H13O_RGR = ({0: ('C', 3), 1: ('C', 3), 2: ('C', 1), 3: ('C', 1), 4: ('C', 1),
               5: ('C', 1), 6: ('C', 2), 7: ('C', 1), 8: ('O', 0)},
              {frozenset({0, 2}): 1, frozenset({1, 3}): 1,
               frozenset({2, 4}): 2, frozenset({3, 5}): 2,
               frozenset({4, 6}): 1, frozenset({5, 7}): 1,
               frozenset({6, 7}): 1, frozenset({8, 7}): 1})
C8H13O_SGR = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
               3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
               6: ('C', 2, None), 7: ('C', 1, False), 8: ('O', 0, None)},
              {frozenset({0, 2}): None, frozenset({1, 3}): None,
               frozenset({2, 4}): False, frozenset({3, 5}): False,
               frozenset({4, 6}): None, frozenset({5, 7}): None,
               frozenset({6, 7}): None, frozenset({8, 7}): None})
C2H2CL2F2_MM_SGR = ({0: ('C', 1, False), 1: ('C', 1, False),
                     2: ('F', 0, None), 3: ('Cl', 0, None),
                     4: ('F', 0, None), 5: ('Cl', 0, None)},
                    {frozenset({0, 1}): None, frozenset({0, 2}): None,
                     frozenset({0, 3}): None, frozenset({1, 4}): None,
                     frozenset({1, 5}): None})
C2H2CL2F2_MP_SGR = ({0: ('C', 1, False), 1: ('C', 1, True),
                     2: ('F', 0, None), 3: ('Cl', 0, None),
                     4: ('F', 0, None), 5: ('Cl', 0, None)},
                    {frozenset({0, 1}): None, frozenset({0, 2}): None,
                     frozenset({0, 3}): None, frozenset({1, 4}): None,
                     frozenset({1, 5}): None})


def test__conn__from_data():
    """ test graph.conn.from_data
    """
    assert C8H13O_CGR == graph.conn.from_data(
        atm_keys=graph.conn.atom_keys(C8H13O_CGR),
        bnd_keys=graph.conn.bond_keys(C8H13O_CGR),
        atm_sym_dct=graph.conn.atom_symbols(C8H13O_CGR),
        atm_hyd_cnt_dct=graph.conn.atom_hydrogen_counts(C8H13O_CGR)
    )


def test__conn__atom_neighbor_keys():
    """ test graph.conn.atom_neighbor_keys
    """
    assert (graph.conn.atom_neighbor_keys(C8H13O_CGR)
            == {0: (2,), 1: (3,), 2: (0, 4), 3: (1, 5), 4: (2, 6), 5: (3, 7),
                6: (4, 7), 7: (5, 6, 8), 8: (7,)})


def test__conn__branch():
    """ test graph.conn.branch
    """
    assert (graph.conn.branch(C8H13O_CGR, 7, 5) ==
            ({1: ('C', 3), 3: ('C', 1), 5: ('C', 1), 7: ('C', 1)},
             {frozenset({1, 3}): None, frozenset({3, 5}): None,
              frozenset({5, 7}): None}))
    assert (graph.conn.branch(C8H13O_CGR, 7, 6) ==
            ({0: ('C', 3), 2: ('C', 1), 4: ('C', 1), 6: ('C', 2), 7: ('C', 1)},
             {frozenset({0, 2}): None, frozenset({2, 4}): None,
              frozenset({4, 6}): None, frozenset({6, 7}): None}))
    assert (graph.conn.branch(C8H13O_CGR, 7, 8) ==
            ({7: ('C', 1), 8: ('O', 0)}, {frozenset({7, 8}): None}))


def test__conn__isomorphism():
    """ test graph.conn.isomorphism
    """
    cgr = C8H13O_CGR
    natms = len(graph.conn.atoms(cgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        cgr_pmt = graph.conn.relabel(cgr, pmt_dct)
        assert graph.conn.isomorphism(cgr, cgr_pmt) == pmt_dct


def test__res__atom_bond_valences():
    """ test graph.res.atom_bond_valences
    """
    assert (graph.res.atom_bond_valences(C8H13O_RGR) ==
            {0: 4, 1: 4, 2: 4, 3: 4, 4: 4, 5: 4, 6: 4, 7: 4, 8: 1})


def test__res__atom_radical_valences():
    """ test graph.res.atom_radical_valences
    """
    assert (graph.res.atom_radical_valences(C8H13O_RGR) ==
            {8: 1})


def test__res__from_data():
    """ test graph.res.from_data
    """
    assert C8H13O_RGR == graph.res.from_data(
        atm_keys=graph.res.atom_keys(C8H13O_RGR),
        bnd_keys=graph.res.bond_keys(C8H13O_RGR),
        atm_sym_dct=graph.res.atom_symbols(C8H13O_RGR),
        atm_hyd_cnt_dct=graph.res.atom_hydrogen_counts(C8H13O_RGR),
        bnd_ord_dct=graph.res.bond_orders(C8H13O_RGR)
    )


def test__res__inchi():
    """ test graph.res.inchi
    """
    rgr = C8H13O_RGR
    natms = len(graph.conn.atoms(rgr))
    for _ in range(10):
        pmt_dct = dict(enumerate(numpy.random.permutation(natms)))
        rgr_pmt = graph.res.relabel(rgr, pmt_dct)
        assert graph.res.inchi(rgr_pmt) == C8H13O_ICH_NO_STEREO
        assert graph.res.atom_inchi_numbers(rgr_pmt) == pmt_dct


def test__stereo__from_data():
    """ test graph.stereo.from_data
    """
    assert C8H13O_SGR == graph.stereo.from_data(
        atm_keys=graph.stereo.atom_keys(C8H13O_SGR),
        bnd_keys=graph.stereo.bond_keys(C8H13O_SGR),
        atm_sym_dct=graph.stereo.atom_symbols(C8H13O_SGR),
        atm_hyd_cnt_dct=graph.stereo.atom_hydrogen_counts(C8H13O_SGR),
        atm_par_dct=graph.stereo.atom_parities(C8H13O_SGR),
        bnd_par_dct=graph.stereo.bond_parities(C8H13O_SGR)
    )


def test__stereo__reflection():
    """ test graph.stereo.reflection
    """
    assert (graph.stereo.reflection(C8H13O_SGR) ==
            ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
              3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
              6: ('C', 2, None), 7: ('C', 1, True), 8: ('O', 0, None)},
             {frozenset({0, 2}): None, frozenset({1, 3}): None,
              frozenset({2, 4}): False, frozenset({3, 5}): False,
              frozenset({4, 6}): None, frozenset({5, 7}): None,
              frozenset({6, 7}): None, frozenset({8, 7}): None}))


def test__stereo__is_chiral():
    """ test graph.stereo.is_chiral
    """
    assert graph.stereo.is_chiral(C2H2CL2F2_MM_SGR) is False
    assert graph.stereo.is_chiral(C2H2CL2F2_MP_SGR) is True


if __name__ == '__main__':
    test__conn__from_data()
    test__conn__atom_neighbor_keys()
    test__conn__branch()
    test__conn__isomorphism()
    test__res__from_data()
    test__res__atom_bond_valences()
    test__res__atom_radical_valences()
    test__res__inchi()
    test__stereo__from_data()
    test__stereo__reflection()
    test__stereo__is_chiral()
