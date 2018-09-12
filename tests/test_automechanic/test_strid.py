""" test automechanic.strid
"""
from automechanic import strid


def test__smiles():
    """ test strid.multiplicity
    """
    assert strid.smiles('C=C1OOC([O])C1=O_m2') == 'C=C1OOC([O])C1=O'


def test__multiplicity():
    """ test strid.multiplicity
    """
    assert strid.multiplicity('C=C1OOC([O])C1=O_m2') == 2


def test__canonical_reaction_identifier():
    """ test strid.canonical_reaction_identifier
    """
    rid = 'O[O]_m2.C=C[C](C(C=C)C)C_m2>>[O][O]_m3.C=CC(C(C=C)C)C_m1'
    can_rid = 'C=CC(C(C=C)C)C_m1.[O][O]_m3>>C=C[C](C(C=C)C)C_m2.O[O]_m2'
    assert strid.canonical_reaction_identifier(rid) == can_rid
    rid = 'O=[C]O_m2>>[H]_m2.O=C=O_m1'
    can_rid = 'O=C=O_m1.[H]_m2>>O=[C]O_m2'
    assert strid.canonical_reaction_identifier(rid) == can_rid


def test__is_spin_balanced():
    """ test strid.is_spin_balanced
    """
    rid1 = '[H][H]_m1.[O]_m3>>[H]_m2.[OH]_m2'
    rid2 = 'OO_m1.[H]_m2>>O[O]_m2.[H][H]_m1'
    rid3 = '[CH3]_m2.[H]_m2>>[CH2]_m1.[H][H]_m1'
    assert strid.is_spin_balanced(rid1) is True
    assert strid.is_spin_balanced(rid2) is True
    assert strid.is_spin_balanced(rid3) is False


def test__is_radical_radical():
    """ test strid.is_radical_radical
    """
    rid1 = '[H][H]_m1.[O]_m3>>[H]_m2.[OH]_m2'
    rid2 = 'OO_m1.[H]_m2>>O[O]_m2.[H][H]_m1'
    rid3 = '[CH3]_m2.[H]_m2>>[CH2]_m1.[H][H]_m1'
    assert strid.is_radical_radical(rid1) is True
    assert strid.is_radical_radical(rid2) is False
    assert strid.is_radical_radical(rid3) is True


if __name__ == '__main__':
    test__is_radical_radical()
