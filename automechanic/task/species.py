""" tasks that operate on CSVs with species information
"""
from . import tablib as tl
from .. import tab
from .. import mol

EXPAND_STEREO__SPC_ID_KEYS = (tl.SPC_ID_SMI_KEY, tl.SPC_ID_ICH_KEY)


def expand_stereo(spc_id_key, spc_csv, spc_csv_out, logger):
    """ expand stereo for species
    """
    assert spc_id_key in EXPAND_STEREO__SPC_ID_KEYS

    schema_keys = (spc_id_key,)
    schema_typs = (tl.SPC_ID_TYP,)
    logger.info("Reading in {:s}".format(spc_csv))
    spc_tbl = tab.read_csv(spc_csv, keys=schema_keys, typs=schema_typs)

    spc_ids = spc_tbl[spc_id_key]

    for spc_smi in spc_ids:
        spc_ich = mol.sm.inchi(spc_smi, force_stereo=True)
        spc_ich_b = mol.ic.inchi_sublayer(spc_ich, key='b')
        spc_ich_t = mol.ic.inchi_sublayer(spc_ich, key='t')
        b_elems = str.split(spc_ich_b, ',') if spc_ich_b else None
        t_elems = str.split(spc_ich_t, ',') if spc_ich_t else None
        if (b_elems and len(b_elems) > 1) or (t_elems and len(t_elems) > 2):
            spc_ich_ = mol.sm.inchi(spc_smi)
            print(spc_smi)
            print(spc_ich_)
            print(spc_ich)
            print()

    spc_ichs = (spc_ids if spc_id_key == tl.SPC_ID_ICH_KEY else
                [mol.sm.inchi(spc_smi) for spc_smi in spc_ids])

    logger.info(len(spc_ichs))

    logger.info(spc_id_key)
    logger.info(spc_csv_out)
