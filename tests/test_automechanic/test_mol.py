""" test the automechanic.mol module
"""
from automechanic import mol

C8H13O_ICHS = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/b5-3-,6-4-/t8-',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/b5-3-,6-4+/t8-',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/b5-3+,6-4-/t8-',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/b5-3+,6-4+/t8-')
C8H13O_ICH_NO_STEREO = 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3'

C5H10O3_ICH = 'InChI=1S/C5H10O3/c1-3-5(7-3)4(2)8-6/h3-6H,1-2H3/t3-,4+,5+/m0/s1'
C5H10O3_ICH_NO_STEREO = 'InChI=1S/C5H10O3/c1-3-5(7-3)4(2)8-6/h3-6H,1-2H3'

C2H2F2_SMI = 'F/C=C/F'
C2H2F2_SMI_NO_STEREO = 'FC=CF'
C2H2F2_ICH = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+'
C2H2F2_ICH_NO_STEREO = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H'
C2H2F2_ICH_STEREO_UNKNOWN = 'InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1?'
C2H2F2_GEO = (('F', (1.584822920001, -0.748486564300, -0.4271224303432)),
              ('C', (0.619219854789, 0.190165523008, -0.2716389279610)),
              ('C', (-0.635730620967, -0.1839138961594, -0.1803643636082)),
              ('F', (-1.602333181611, 0.736677675476, -0.02605091648865)),
              ('H', (0.916321356258, 1.229945559249, -0.2271265738271)),
              ('H', (-0.882300328471, -1.224388297273, -0.229635969682)))
C2H2F2_ICH_ORDER = (1, 2, 0, 3)

C2H4CLF_SMI = 'C(Cl)(F)C'


def test__ge__inchi_with_order():
    """ test mol.ge.inchi_with_order
    """
    assert (mol.ge.inchi_with_order(C2H2F2_GEO)
            == (C2H2F2_ICH, C2H2F2_ICH_ORDER))


def test__ge__inchi():
    """ test mol.ge.inchi
    """
    assert mol.ge.inchi(C2H2F2_GEO) == C2H2F2_ICH


def test__sm__inchi():
    """ test mol.sm.inchi
    """
    assert mol.sm.inchi(C2H2F2_SMI) == C2H2F2_ICH
    assert mol.sm.inchi(C2H2F2_SMI_NO_STEREO) == C2H2F2_ICH_NO_STEREO
    assert (mol.sm.inchi(C2H2F2_SMI_NO_STEREO, force_stereo=True)
            == C2H2F2_ICH_STEREO_UNKNOWN)

    assert mol.sm.inchi(C2H4CLF_SMI) == 'InChI=1S/C2H4ClF/c1-2(3)4/h2H,1H3'
    assert (mol.sm.inchi(C2H4CLF_SMI, force_stereo=True)
            == 'InChI=1/C2H4ClF/c1-2(3)4/h2H,1H3/t2?')


def test__ic__inchi():
    """ test mol.ic.inchi
    """
    assert mol.ic.inchi(C2H2F2_ICH_NO_STEREO) == C2H2F2_ICH_NO_STEREO
    assert (mol.ic.inchi(C2H2F2_ICH_NO_STEREO, force_stereo=True)
            == C2H2F2_ICH_STEREO_UNKNOWN)


def test__ic__inchi_prefix():
    """ test mol.ic.inchi_prefix
    """
    assert mol.ic.inchi_prefix(C2H2F2_ICH) == 'InChI=1S'
    assert mol.ic.inchi_prefix(C2H2F2_ICH_NO_STEREO) == 'InChI=1S'
    assert mol.ic.inchi_prefix(C2H2F2_ICH_STEREO_UNKNOWN) == 'InChI=1'


def test__ic__inchi_formula():
    """ test mol.ic.inchi_formula
    """
    assert mol.ic.inchi_formula(C2H2F2_ICH) == 'C2H2F2'


def test__ic__inchi_sublayer():
    """ test mol.ic.inchi_sublayer
    """
    assert mol.ic.inchi_sublayer(C2H2F2_ICH, 'c') == '3-1-2-4'
    assert mol.ic.inchi_sublayer(C2H2F2_ICH, 'h') == '1-2H'
    assert mol.ic.inchi_sublayer(C2H2F2_ICH, 'b') == '2-1+'
    assert mol.ic.inchi_sublayer(C2H2F2_ICH_STEREO_UNKNOWN, 'c') == '3-1-2-4'
    assert mol.ic.inchi_sublayer(C2H2F2_ICH_STEREO_UNKNOWN, 'h') == '1-2H'
    assert mol.ic.inchi_sublayer(C2H2F2_ICH_STEREO_UNKNOWN, 'b') == '2-1?'
    assert mol.ic.inchi_sublayer(C2H2F2_ICH_NO_STEREO, 'c') == '3-1-2-4'
    assert mol.ic.inchi_sublayer(C2H2F2_ICH_NO_STEREO, 'h') == '1-2H'
    assert mol.ic.inchi_sublayer(C2H2F2_ICH_NO_STEREO, 'b') is None


def test__ic__normalized_inchi():
    """ test mol.ic.normalized_inchi
    """
    assert (mol.ic.normalized_inchi(C2H2F2_ICH)
            == 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+')
    assert (mol.ic.normalized_inchi(C2H2F2_ICH_NO_STEREO)
            == 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H')
    assert (mol.ic.normalized_inchi(C2H2F2_ICH_STEREO_UNKNOWN)
            == 'InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1?')
    assert (mol.ic.normalized_inchi(C8H13O_ICH_NO_STEREO)
            == 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3')
    assert (mol.ic.normalized_inchi(C8H13O_ICHS[0])
            == 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
               'b5-3-,6-4-/t8-')


def test__ic__inchi_expand_unknown_stereo():
    """ test mol.ic.inchi_expand_unknown_stereo
    """
    assert (mol.ic.inchi_expand_unknown_stereo(C8H13O_ICH_NO_STEREO)
            == C8H13O_ICHS)
    assert (mol.ic.inchi_expand_unknown_stereo(C8H13O_ICHS[0])
            == (C8H13O_ICHS[0],))


def test__ic__inchi_key():
    """ test mol.ic.inchi_key
    """
    assert mol.ic.inchi_key(C2H2F2_ICH) == 'WFLOTYSKFUPZQB-OWOJBTEDSA-N'
    assert (mol.ic.inchi_key(C2H2F2_ICH_NO_STEREO)
            == 'WFLOTYSKFUPZQB-UHFFFAOYSA-N')


def test__ic__geometry():
    """ test mol.ic.geometry
    """
    # make sure these run
    mol.ic.geometry(C2H2F2_ICH)
    for ich in C8H13O_ICHS:
        mol.ic.geometry(ich)


def test__ge__atoms():
    """ test mol.ge.atoms()
    """
    assert mol.ge.atoms(C2H2F2_GEO) == ('F', 'C', 'C', 'F', 'H', 'H')


def test__ge__bonds():
    """ test mol.ge.bonds()
    """
    assert (mol.ge.bonds(C2H2F2_GEO)
            == {frozenset({0, 1}): 1, frozenset({1, 2}): 2,
                frozenset({1, 4}): 1, frozenset({2, 3}): 1,
                frozenset({2, 5}): 1})


def test__ge__graph():
    """ test mol.ge.graph()
    """
    assert (mol.ge.graph(C2H2F2_GEO)
            == (('F', 'C', 'C', 'F', 'H', 'H'),
                {frozenset({0, 1}): 1, frozenset({1, 2}): 2,
                 frozenset({1, 4}): 1, frozenset({2, 3}): 1,
                 frozenset({2, 5}): 1}))


def test__gr__atoms():
    """ test mol.gr.atoms()
    """
    gra = (('C', 'C', 'H', 'H', 'H', 'H'),
           {frozenset([0, 3]): 1, frozenset([1, 4]): 1, frozenset([0, 2]): 1,
            frozenset([1, 5]): 1, frozenset([0, 1]): 2})
    assert mol.gr.atoms(gra) == ('C', 'C', 'H', 'H', 'H', 'H')


def test__gr__bonds():
    """ test mol.gr.bonds()
    """
    gra = (('C', 'C', 'H', 'H', 'H', 'H'),
           {frozenset([0, 3]): 1, frozenset([1, 4]): 1, frozenset([0, 2]): 1,
            frozenset([1, 5]): 1, frozenset([0, 1]): 2})
    assert (mol.gr.bonds(gra) ==
            {frozenset([0, 3]): 1, frozenset([1, 4]): 1, frozenset([0, 2]): 1,
             frozenset([1, 5]): 1, frozenset([0, 1]): 2})


def test__gr__indices():
    """ test mol.gr.indices()
    """
    gra = (('C', 'C', 'H', 'H', 'H', 'H'),
           {frozenset([0, 3]): 1, frozenset([1, 4]): 1, frozenset([0, 2]): 1,
            frozenset([1, 5]): 1, frozenset([0, 1]): 2})
    assert mol.gr.indices(gra) == (0, 1, 2, 3, 4, 5)


def test__gr__atom_at():
    """ test mol.gr.atom_at()
    """
    gra = (('C', 'C', 'H', 'H', 'H', 'H'),
           {frozenset([0, 3]): 1, frozenset([1, 4]): 1, frozenset([0, 2]): 1,
            frozenset([1, 5]): 1, frozenset([0, 1]): 2})
    assert mol.gr.atom_at(gra, 0) == 'C'
    assert mol.gr.atom_at(gra, 1) == 'C'
    assert mol.gr.atom_at(gra, 2) == 'H'


def test__gr__bonds_at():
    """ test mol.gr.bonds_at()
    """
    gra = (('C', 'C', 'H', 'H', 'H', 'H'),
           {frozenset([0, 3]): 1, frozenset([1, 4]): 1, frozenset([0, 2]): 1,
            frozenset([1, 5]): 1, frozenset([0, 1]): 2})
    assert (mol.gr.bonds_at(gra, 0) ==
            {frozenset({0, 3}): 1, frozenset({0, 2}): 1, frozenset({0, 1}): 2})


def test__gr__bound_electrons_at():
    """ test mol.gr.bound_electrons_at()
    """
    gra = (('C', 'C', 'H', 'H', 'H', 'H'),
           {frozenset([0, 3]): 1, frozenset([1, 4]): 1, frozenset([0, 2]): 1,
            frozenset([1, 5]): 1, frozenset([0, 1]): 2})
    assert mol.gr.bound_electrons_at(gra, 0) == 4
    assert mol.gr.bound_electrons_at(gra, 1) == 4
    assert mol.gr.bound_electrons_at(gra, 2) == 1


def test__gr__radical_electrons_at():
    """ test mol.gr.radical_electrons_at()
    """
    gra = (('C', 'C', 'H', 'H'),
           {frozenset([0, 3]): 1, frozenset([0, 2]): 1, frozenset([0, 1]): 2})
    assert mol.gr.radical_electrons_at(gra, 0) == 0
    assert mol.gr.radical_electrons_at(gra, 1) == 2
    assert mol.gr.radical_electrons_at(gra, 2) == 0
    assert mol.gr.radical_electrons_at(gra, 3) == 0


def test__gr__radicals():
    """ test mol.gr.radicals()
    """
    gra = (('C', 'C', 'H', 'H'),
           {frozenset([0, 3]): 1, frozenset([0, 2]): 1, frozenset([0, 1]): 2})
    assert mol.gr.radicals(gra) == {1: 2}


if __name__ == '__main__':
    test__ge__inchi_with_order()
    test__sm__inchi()
    test__ic__inchi()
    test__ic__inchi_prefix()
    test__ic__inchi_formula()
    test__ic__inchi_sublayer()
    test__ic__normalized_inchi()
    test__ic__inchi_expand_unknown_stereo()
    test__ic__inchi_key()
    test__ic__geometry()
    test__ge__atoms()
    test__ge__bonds()
    test__ge__graph()
    test__gr__atoms()
    test__gr__bonds()
    test__gr__indices()
    test__gr__atom_at()
    test__gr__bonds_at()
    test__gr__bound_electrons_at()
    test__gr__radical_electrons_at()
    test__gr__radicals()
