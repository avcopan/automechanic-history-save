""" test the automechanic.pchemkin module
"""
import os
from automechanic import pchemkin2 as pchemkin

PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')


def test__reactions():
    """ test pchemkin.reactions
    """
    mech_fpath = os.path.join(DATA_PATH, 'heptane_mechanism.txt')
    mech_str = open(mech_fpath).read()
    reacs = pchemkin.reactions(mech_str)
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


if __name__ == '__main__':
    test__reactions()
