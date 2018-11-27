""" tasks that operate on CSVs with species information
"""
from . import tablib as tl
from .. import tab

GEOM_SPEC_KEYS = (
    tl.GEOM_SPEC_SMI_KEY,
    # tl.GEOM_SPEC_ICH_KEY,
    # tl.GEOM_SPEC_XYZ_KEY,
)


def fill_guess_geoms(geom_type_key, spc_csv, spc_csv_out, db_prefix, logger):
    """ determine guess geometries for structure optimization
    """
    assert geom_type_key in GEOM_SPEC_KEYS

    schema_keys = (geom_type_key, tl.SPC_MULT_KEY)
    schema_typs = (tl.GEOM_SPEC_TYP, tl.SPC_MULT_TYP)

    logger.info("Reading in {:s}".format(spc_csv))
    spc_tbl = tab.read_csv(spc_csv, keys=schema_keys, typs=schema_typs)

    geom_specs = spc_tbl[geom_type_key]
    geom_mults = spc_tbl[tl.SPC_MULT_KEY]

    logger.info(list(geom_specs))
    logger.info(list(geom_mults))

    logger.info(geom_type_key)
    logger.info(spc_csv_out)
    logger.info(db_prefix)
