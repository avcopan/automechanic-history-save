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
