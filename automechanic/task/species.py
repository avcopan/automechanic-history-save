""" tasks that operate on CSVs with species information
"""
from .. import params as par
from .. import tab
from .. import mol
from .. import fslib
from .. import fs
from ..iohelp import timestamp_if_exists


TO_INCHI__SPC_ID_VALS = (
    par.SPC_ID_SMI_KEY, par.SPC_ID_ICH_KEY)

FILESYSTEM__STEREO_HANDLING_VALS = (
    par.EXPAND_STEREO, par.PICK_STEREO)

STEREO_HANDLING_DEF = par.PICK_STEREO


def to_inchi(spc_id_key, spc_csv, spc_csv_out, logger):
    """ convert species identifiers to InChI
    """
    assert spc_id_key in TO_INCHI__SPC_ID_VALS

    logger.info("Reading in {:s}".format(spc_csv))
    tbl = tab.read_csv(spc_csv)

    logger.info("Converting to InChI")
    tbl = _to_inchi(tbl, spc_id_key)

    logger.info("Writing to {:s}".format(spc_csv_out))
    timestamp_if_exists(spc_csv_out)
    tab.write_csv(spc_csv_out, tbl)


def filesystem(spc_csv, spc_csv_out, stereo_handling, filesystem_prefix,
               logger):
    """ chart the species filesystem structure
    """
    assert stereo_handling in FILESYSTEM__STEREO_HANDLING_VALS

    logger.info("Reading in {:s}".format(spc_csv))
    tbl = tab.read_csv(spc_csv)

    logger.info("Handling stereo in mode '{:s}'".format(stereo_handling))
    tbl = _handle_stereo(tbl, mode=stereo_handling)

    logger.info("Creating filesystem at '{:s}'".format(filesystem_prefix))
    tbl = _create_filesystem(tbl, fs_root_pth=filesystem_prefix)

    logger.info("Writing to {:s}".format(spc_csv_out))
    timestamp_if_exists(spc_csv_out)
    tab.write_csv(spc_csv_out, tbl)

    logger.info(filesystem_prefix)


def _to_inchi(tbl, spc_id_key):
    assert spc_id_key in (par.SPC_ID_SMI_KEY, par.SPC_ID_ICH_KEY)
    tbl = tab.enforce_schema(tbl, keys=(spc_id_key,), typs=(par.SPC_ID_TYP,))

    conv_ = (mol.inchi.recalculate if spc_id_key == par.SPC_ID_ICH_KEY else
             mol.smiles.inchi)

    sids = tbl[spc_id_key]
    ichs = list(map(conv_, sids))
    tbl = tbl[[key for key in tab.keys_(tbl) if key != spc_id_key]]
    tbl[par.SPC_ID_ICH_KEY] = ichs
    return tbl


def _handle_stereo(tbl, mode):
    assert mode in (par.EXPAND_STEREO, par.PICK_STEREO)
    return (_handle_stereo_by_expanding(tbl) if mode == par.EXPAND_STEREO else
            _handle_stereo_by_picking(tbl))


def _handle_stereo_by_picking(tbl):
    tbl = tab.enforce_schema(tbl,
                             keys=(par.SPC_ID_ICH_KEY,),
                             typs=(par.SPC_ID_TYP,))
    tbl = tbl.copy()
    # use coordinates to get stereo assignments
    ichs = tbl[par.SPC_ID_ICH_KEY]
    ichs = list(map(mol.geom.inchi, map(mol.inchi.geometry, ichs)))
    assert not any(map(mol.inchi.has_unknown_stereo_elements, ichs))
    tbl[par.SPC_ID_ICH_KEY] = ichs
    return tbl


def _handle_stereo_by_expanding(tbl):
    tbl = tab.enforce_schema(tbl,
                             keys=(par.SPC_ID_ICH_KEY,),
                             typs=(par.SPC_ID_TYP,))

    # do the stereo expansion
    ichs = tbl[par.SPC_ID_ICH_KEY]
    ichsts_lst = list(map(mol.inchi.compatible_stereoisomers, ichs))

    # save the current index before expansion
    idx_save_key = tab.next_index_save_key(tbl)
    tbl = tab.save_index(tbl)

    # make the stereo expansion table
    tbl_idxs = tab.idxs_(tbl)
    vals = [[idx, ichst]
            for idx, ichsts in zip(tbl_idxs, ichsts_lst) for ichst in ichsts]
    keys = (idx_save_key, par.SPC_ID_ICH_KEY)
    typs = (tab.IDX_TYP, par.SPC_ID_TYP)
    ste_tbl = tab.from_records(vals, keys=keys, typs=typs)

    # make the return table
    keys = [key for key in tab.keys_(tbl) if key != par.SPC_ID_ICH_KEY]
    tbl = tab.left_join(tbl[keys], ste_tbl, key=idx_save_key)
    return tbl


def _create_filesystem(tbl, fs_root_pth):
    id_keys = (par.SPC_ID_ICH_KEY, par.SPC_MULT_KEY)
    id_typs = (par.SPC_ID_TYP, par.SPC_MULT_TYP)
    tbl = tab.enforce_schema(tbl, keys=id_keys, typs=id_typs)

    def __create_branch(ich, mult):
        print((ich, mult))
        sgms = fslib.species.branch_segments(ich, mult)
        return fs.branch.create(sgms)

    with fs.enter(fs_root_pth):
        pth_tbl = tab.from_starmap(tbl, __create_branch, id_keys,
                                   keys=(par.SPC_FS_PATH_KEY,),
                                   typs=(par.SPC_FS_PATH_TYP,))
        tbl = tab.left_join(tbl, pth_tbl)

    return tbl
