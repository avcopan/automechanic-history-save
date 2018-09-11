""" test the automechanic.pchemkin module
"""
import os
from automechanic import pchemkin

PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')


def test__reactions():
    """ test pchemkin.reactions
    """
    mech_fpath = os.path.join(DATA_PATH, 'heptane_mechanism.txt')
    mech_str = open(mech_fpath).read()
    reacs = pchemkin.reactions_without_em(mech_str)
    assert len(reacs) == 4


def test__remove_comments():
    """ test pchemkin.remove_comments
    """
    mech_fpath = os.path.join(DATA_PATH, 'mechanism.txt')
    clean_mech_fpath = os.path.join(DATA_PATH, 'clean_mechanism.txt')
    mech_str = open(mech_fpath).read()
    clean_mech_str = open(clean_mech_fpath).read()
    assert pchemkin.remove_comments(mech_str) == clean_mech_str


def test__block():
    """ test pchemkin.block
    """
    mech_fpath = os.path.join(DATA_PATH, 'mechanism.txt')
    mech_str = open(mech_fpath).read()

    elem_block_fpath = os.path.join(DATA_PATH, 'elem_block.txt')
    spec_block_fpath = os.path.join(DATA_PATH, 'spec_block.txt')
    ther_block_fpath = os.path.join(DATA_PATH, 'ther_block.txt')
    reac_block_fpath = os.path.join(DATA_PATH, 'reac_block.txt')
    elem_block_str = open(elem_block_fpath).read()
    spec_block_str = open(spec_block_fpath).read()
    ther_block_str = open(ther_block_fpath).read()
    reac_block_str = open(reac_block_fpath).read()
    assert pchemkin.elements_block(mech_str).strip() == elem_block_str.strip()
    assert pchemkin.species_block(mech_str).strip() == spec_block_str.strip()
    assert pchemkin.thermo_block(mech_str).strip() == ther_block_str.strip()
    assert pchemkin.reactions_block(mech_str).strip() == reac_block_str.strip()


def test__species():
    """ test pchemkin.species
    """
    mech_fpath = os.path.join(DATA_PATH, 'mechanism.txt')
    mech_str = open(mech_fpath).read()
    specs = ('N2', 'Ne', 'CO(1)', 'H2(2)', 'O2(3)', 'H(4)', 'O(5)', 'OH(6)',
             'H2O(7)', 'HO2(10)', 'H2O2(11)', 'CO2(12)', 'HCO(13)',
             'CHO2(38)', 'HOCO(39)', 'CHO4(44)', 'CHO4(45)', 'CHO3(61)',
             'CH2O3(82)')
    assert pchemkin.species(mech_str) == specs


def test__reactions_without_em():
    """ test pchemkin.reactions_without_em
    """
    mech_fpath = os.path.join(DATA_PATH, 'mechanism.txt')
    mech_str = open(mech_fpath).read()
    reacs = ((('H(4)', 'HO2(10)'), ('H2(2)', 'O2(3)')),
             (('CO(1)', 'O2(3)'), ('O(5)', 'CO2(12)')),
             (('O2(3)', 'H(4)'), ('O(5)', 'OH(6)')),
             (('H2(2)', 'O(5)'), ('H(4)', 'OH(6)')),
             (('H2(2)', 'O(5)'), ('H(4)', 'OH(6)')),
             (('H(4)', 'HO2(10)'), ('OH(6)', 'OH(6)')),
             (('O(5)', 'HO2(10)'), ('O2(3)', 'OH(6)')),
             (('CO(1)', 'OH(6)'), ('H(4)', 'CO2(12)')),
             (('CO(1)', 'OH(6)'), ('H(4)', 'CO2(12)')),
             (('CO(1)', 'HO2(10)'), ('OH(6)', 'CO2(12)')),
             (('H2(2)', 'OH(6)'), ('H(4)', 'H2O(7)')),
             (('OH(6)', 'OH(6)'), ('O(5)', 'H2O(7)')),
             (('H2O(7)', 'H2O(7)'), ('H(4)', 'OH(6)', 'H2O(7)')),
             (('H(4)', 'HO2(10)'), ('O(5)', 'H2O(7)')),
             (('OH(6)', 'HO2(10)'), ('O2(3)', 'H2O(7)')),
             (('OH(6)', 'HO2(10)'), ('O2(3)', 'H2O(7)')),
             (('H(4)', 'HCO(13)'), ('CO(1)', 'H2(2)')),
             (('O(5)', 'HCO(13)'), ('CO(1)', 'OH(6)')),
             (('O(5)', 'HCO(13)'), ('H(4)', 'CO2(12)')),
             (('OH(6)', 'HCO(13)'), ('CO(1)', 'H2O(7)')),
             (('O2(3)', 'HCO(13)'), ('CO(1)', 'HO2(10)')),
             (('HO2(10)', 'HO2(10)'), ('O2(3)', 'H2O2(11)')),
             (('HO2(10)', 'HO2(10)'), ('O2(3)', 'H2O2(11)')),
             (('H(4)', 'H2O2(11)'), ('OH(6)', 'H2O(7)')),
             (('H(4)', 'H2O2(11)'), ('H2(2)', 'HO2(10)')),
             (('O(5)', 'H2O2(11)'), ('OH(6)', 'HO2(10)')),
             (('OH(6)', 'H2O2(11)'), ('H2O(7)', 'HO2(10)')),
             (('OH(6)', 'H2O2(11)'), ('H2O(7)', 'HO2(10)')),
             (('H(4)', 'HOCO(39)'), ('H2(2)', 'CO2(12)')),
             (('O2(3)', 'HOCO(39)'), ('HO2(10)', 'CO2(12)')),
             (('OH(6)', 'HOCO(39)'), ('H2O(7)', 'CO2(12)')),
             (('O(5)', 'HOCO(39)'), ('OH(6)', 'CO2(12)')),
             (('HO2(10)', 'HOCO(39)'), ('H2O2(11)', 'CO2(12)')),
             (('H(4)', 'CHO2(38)'), ('H2(2)', 'CO2(12)')),
             (('O2(3)', 'CHO2(38)'), ('HO2(10)', 'CO2(12)')),
             (('OH(6)', 'CHO2(38)'), ('H2O(7)', 'CO2(12)')),
             (('O(5)', 'CHO2(38)'), ('OH(6)', 'CO2(12)')),
             (('HO2(10)', 'CHO2(38)'), ('H2O2(11)', 'CO2(12)')),
             (('OH(6)', 'CH2O3(82)'), ('H2O(7)', 'CHO3(61)')),
             (('HO2(10)', 'CHO3(61)'), ('O2(3)', 'CH2O3(82)')),
             (('H2O2(11)', 'CHO3(61)'), ('HO2(10)', 'CH2O3(82)')),
             (('CHO2(38)', 'CHO3(61)'), ('CO2(12)', 'CH2O3(82)')),
             (('HOCO(39)', 'CHO3(61)'), ('CO2(12)', 'CH2O3(82)')))
    assert pchemkin.reactions_without_em(mech_str) == reacs


def test__reactions_with_em():
    """ test pchemkin.reactions_with_em
    """
    mech_fpath = os.path.join(DATA_PATH, 'mechanism.txt')
    mech_str = open(mech_fpath).read()
    reacs = ((('H2(2)',), ('H(4)', 'H(4)')),
             (('O(5)', 'O(5)'), ('O2(3)',)),
             (('H(4)', 'O(5)'), ('OH(6)',)),
             (('H2O(7)',), ('H(4)', 'OH(6)')),
             (('HCO(13)',), ('CO(1)', 'H(4)')))
    assert pchemkin.reactions_with_em(mech_str) == reacs


def test__reactions_with_parentheses_em():
    """ test pchemkin.reactions_with_parentheses_em
    """
    mech_fpath = os.path.join(DATA_PATH, 'mechanism.txt')
    mech_str = open(mech_fpath).read()
    reacs = ((('O2(3)', 'H(4)'), ('HO2(10)',)),
             (('CO(1)', 'O(5)'), ('CO2(12)',)),
             (('O(5)', 'OH(6)'), ('HO2(10)',)),
             (('H2O2(11)',), ('OH(6)', 'OH(6)')),
             (('H(4)', 'HO2(10)'), ('H2O2(11)',)),
             (('H(4)', 'CO2(12)'), ('O(5)', 'HCO(13)')),
             (('CO(1)', 'HO2(10)'), ('OH(6)', 'CO2(12)')),
             (('CHO4(45)',), ('HO2(10)', 'CO2(12)')),
             (('HO2(10)', 'CO2(12)'), ('O2(3)', 'HOCO(39)')),
             (('HOCO(39)',), ('H(4)', 'CO2(12)')),
             (('OH(6)', 'CO2(12)'), ('O(5)', 'HOCO(39)')),
             (('CO(1)', 'OH(6)'), ('HOCO(39)',)),
             (('CHO4(45)',), ('O2(3)', 'HOCO(39)')),
             (('HO2(10)', 'CO2(12)'), ('O(5)', 'CHO3(61)')),
             (('CHO3(61)',), ('OH(6)', 'CO2(12)')),
             (('CO(1)', 'HO2(10)'), ('CHO3(61)',)),
             (('CHO4(45)',), ('O(5)', 'CHO3(61)')),
             (('CHO3(61)',), ('O(5)', 'HOCO(39)')),
             (('CHO2(38)',), ('H(4)', 'CO2(12)')),
             (('O2(3)', 'HCO(13)'), ('O(5)', 'CHO2(38)')),
             (('CHO2(38)',), ('HOCO(39)',)),
             (('CHO2(38)',), ('O(5)', 'HCO(13)')),
             (('CHO4(44)',), ('HO2(10)', 'CO2(12)')),
             (('CHO4(44)',), ('CHO4(45)',)),
             (('CO(1)', 'H2O(7)'), ('OH(6)', 'HCO(13)')),
             (('CO(1)', 'H2O(7)'), ('H(4)', 'HOCO(39)')),
             (('CO(1)', 'H2O(7)'), ('H(4)', 'CHO2(38)')),
             (('CO(1)', 'H2O(7)'), ('H2(2)', 'CO2(12)')),
             (('CH2O3(82)',), ('OH(6)', 'HOCO(39)')),
             (('CH2O3(82)',), ('H(4)', 'CHO3(61)')))
    assert pchemkin.reactions_with_parentheses_em(mech_str) == reacs


def test__thermo_functions():
    """ test thermo functions
    """
    import numpy
    from automechanic import therm

    mech_fpath = os.path.join(DATA_PATH, 'mechanism.txt')
    mech_str = open(mech_fpath).read()
    cfts_lo_dct = pchemkin.thermo_lower_coefficients(mech_str)
    cfts_hi_dct = pchemkin.thermo_upper_coefficients(mech_str)
    tcom_dct = pchemkin.thermo_common_temperatures(mech_str)
    for spc, tcom in tcom_dct.items():
        cfts_lo = cfts_lo_dct[spc]
        cfts_hi = cfts_hi_dct[spc]
        tlo = numpy.nextafter(tcom, tcom-1.)
        thi = numpy.nextafter(tcom, tcom+1.)
        enth_lo = therm.enthalpy(tlo, cfts_lo)
        enth_hi = therm.enthalpy(thi, cfts_hi)
        assert numpy.allclose(enth_lo, enth_hi)


if __name__ == '__main__':
    test__reactions()
