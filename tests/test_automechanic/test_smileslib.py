""" test the automechanic.smileslib module
"""
from automechanic import smileslib


def test__match_hydrogen_abstraction_formula():
    """ test smileslib.match_hydrogen_abstraction_formula()
    """
    smi = '[CH3].O>>[OH].C'
    ref_smi = 'C.[OH]>>[CH3].O'
    assert smileslib.match_hydrogen_abstraction_formula(smi) == ref_smi
    smi = '[CH2]C=CC=C.[H]>>CC=CC=C'
    assert smileslib.match_hydrogen_abstraction_formula(smi) is None
    smi = 'CC=C.C=C=C>>[CH2]C=C.[CH2]C=C'
    ref_smi = 'CC=C.C=C=C>>[CH2]C=C.[CH2]C=C'
    assert smileslib.match_hydrogen_abstraction_formula(smi) == ref_smi


if __name__ == '__main__':
    test__match_hydrogen_abstraction_formula()
