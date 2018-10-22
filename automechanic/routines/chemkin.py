""" CHEMKIN routines
"""
# from .io import read_txt
# from ..parse.chemkin import species_names
# from ..parse.chemkin import thermo_data
# from ..parse.chemkin import reaction_data

# species keys
SPC_KEY = 'species'
NASA_LO_KEYS = ('nasa_lo_1', 'nasa_lo_2', 'nasa_lo_3', 'nasa_lo_4',
                'nasa_lo_5', 'nasa_lo_6', 'nasa_lo_7')
NASA_HI_KEYS = ('nasa_hi_1', 'nasa_hi_2', 'nasa_hi_3', 'nasa_hi_4',
                'nasa_hi_5', 'nasa_hi_6', 'nasa_hi_7')
NASA_TC_KEY = 'nasa_tc'

# reaction keys
PRDS_KEY = 'products'
RCTS_KEY = 'reactants'
ARRH_KEYS = ('arrh_a', 'arrh_b', 'arrh_b')


def to_csv(mech_txt_lst, rxn_csv_out, spc_csv_out, logger):
    """ parse CHEMKIN information to CSV
    """
    logger.info("Reading in mechanism file(s)")
    # mech_str = str.join('\n', map(read_txt, mech_txt_lst))

    logger.info("Finding species data")
    # _species_table(mech_str)

    logger.info("Finding reactions data")
    # _reactions_table(mech_str)

    logger.debug(mech_txt_lst)
    logger.debug(rxn_csv_out)
    logger.debug(spc_csv_out)


# def _species_table(mech_str):
#     spcs = species_names(mech_str)
#     thm_dat_lst = thermo_data(mech_str)
#     assert len(thm_dat_lst) == len(spcs)
#     spcs, nasa_los, nasa_his, nasa_tcs = zip(*thm_dat_lst)


# def _reactions_table(mech_str):
#     rct_dat_lst = reaction_data(mech_str)
#     print(rct_dat_lst)
