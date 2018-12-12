""" tasks that operate on CHEMKIN-format files
"""
from .. import params as par
from .. import tab
from ..iohelp import read_string
from ..iohelp import timestamp_if_exists
from ..parse.chemkin import species_names
from ..parse.chemkin import thermo_data
from ..parse.chemkin import reaction_data


def to_csv(mech_txt_lst, rxn_csv_out, spc_csv_out, logger):
    """ parse CHEMKIN information to CSV
    """
    logger.info("Reading in mechanism file(s)")
    mech_str = '\n'.join(map(read_string, mech_txt_lst))

    logger.info("Finding species data")
    spc_tbl = _species_table(mech_str)

    logger.info("Writing species data to {:s}".format(spc_csv_out))
    timestamp_if_exists(spc_csv_out)
    tab.write_csv(spc_csv_out, spc_tbl, float_format='%.8f')

    logger.info("Finding reactions data")
    rxn_tbl = _reactions_table(mech_str)

    logger.info("Writing reaction data to {:s}".format(spc_csv_out))
    timestamp_if_exists(rxn_csv_out)
    tab.write_csv(rxn_csv_out, rxn_tbl, float_format='%.8f')


def _species_table(mech_str):
    spcs = species_names(mech_str)
    thm_dat_lst = thermo_data(mech_str)
    assert len(thm_dat_lst) == len(spcs)
    typs = (par.SPC_NAME_TYP, par.NASA_C_TYP, par.NASA_C_TYP, par.NASA_T_TYP)
    keys = (par.SPC_NAME_KEY, par.NASA_C_LO_KEYS, par.NASA_C_HI_KEYS,
            par.NASA_T_KEYS)
    spc_tbl = tab.from_records(vals=thm_dat_lst, keys=keys, typs=typs)
    return spc_tbl


def _reactions_table(mech_str):
    rxn_dat_lst = reaction_data(mech_str)
    typs = (par.RXN_NAME_TYP, par.ARRH_TYP)
    keys = (par.RXN_NAME_KEY, par.ARRH_KEYS)
    rxn_tbl = tab.from_records(vals=rxn_dat_lst, keys=keys, typs=typs)
    return rxn_tbl
