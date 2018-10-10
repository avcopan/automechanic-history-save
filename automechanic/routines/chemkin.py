""" CHEMKIN routines
"""
from .io import read_txt
from ..parse.chemkin import species_names
from ..parse.chemkin import thermo_data
from ..parse.chemkin import reaction_data


def to_csv(mech_txt_lst, rxn_csv_out, spc_csv_out, logger):
    """ parse CHEMKIN information to CSV
    """
    logger.info("Reading in mechanism file(s)")
    mech_str = str.join('\n', map(read_txt, mech_txt_lst))
    spcs = species_names(mech_str)
    thm_dat_lst = thermo_data(mech_str)
    rct_dat_lst = reaction_data(mech_str)

    logger.debug(len(spcs))
    logger.debug(len(thm_dat_lst))
    logger.debug(len(rct_dat_lst))
    logger.debug(rxn_csv_out)
    logger.debug(spc_csv_out)
