""" drivers and functions for I/O
"""
import os
import time
import json
import subprocess
from functools import partial
import pandas
from .iohelp2 import addition_indices
from .iohelp2 import addition_candidate
from .iohelp2 import addition_xyz_string_dictionary
from .iohelp2 import addition_input_string
from .iohelp2 import abstraction_indices
from .iohelp2 import abstraction_candidate
from .iohelp2 import abstraction_xyz_string_dictionary
from .iohelp2 import abstraction_input_string


def timestamp_if_exists(fpath):
    """ open a file, avoiding overwrites if requested
    """
    if os.path.exists(fpath):
        time_stamp = time.strftime("%Y%m%d-%H%M%S")
        fpath = "{:s}_{:s}".format(fpath, time_stamp)
    return fpath


def write_file(fpath, contents, mode='w'):
    """ write contents to a file
    """
    fle = open(fpath, mode)
    fle.write(contents)
    fle.close()


def read_file(fpath):
    """ read file contents as a string
    """
    return open(fpath).read()


def init_geometries(spc_csv, spc_csv_out, geom_dir, sid2fname,
                    logger):
    """ initialize species geometries
    """
    from .strid import xyz_string

    logger.info("Initializing species geometries")

    spc_df = pandas.read_csv(spc_csv)

    assert 'species_id' in spc_df

    if not os.path.exists(geom_dir):
        os.mkdir(geom_dir)

    for idx, row in spc_df.iterrows():
        sid = row['species_id']
        fname = sid2fname(sid)
        fpath = os.path.join(geom_dir, fname)
        if not os.path.exists(fpath):
            logger.debug("Writing geometry for {:s} to {:s}"
                         .format(sid, fpath))
            dxyz = xyz_string(sid)
            write_file(fpath, contents=dxyz)
        else:
            logger.debug("Geometry for {:s} already exists at {:s}"
                         .format(sid, fpath))
        spc_df.at[idx, 'path'] = fpath

    logger.info("Writing species .xyz paths to {:s}".format(spc_csv_out))
    spc_df.to_csv(spc_csv_out, index=False)


def geometries(spc_csv):
    """ a dictionary of geometries, indexed by species ID
    """
    from .dotxyz import geometry

    prefix = os.path.dirname(spc_csv)
    add_prefix_ = partial(os.path.join, prefix)
    spc_df = pandas.read_csv(spc_csv)
    sids = spc_df['species_id']
    paths = map(add_prefix_, spc_df['path'])
    dxyzs = map(read_file, paths)
    mgeos = map(geometry, dxyzs)
    mgeo_dct = dict(zip(sids, mgeos))
    return mgeo_dct


def init_from_rmg(rmg_mech_json, spc_csv_out, rxn_csv_out, logger):
    """ initialize the mechanism from RMG's JSON file
    """
    from .prmg import mechanism_species_identifiers
    from .prmg import mechanism_reaction_identifiers
    from .prmg import mechanism_uncertainties
    from .prmg import mechanism_sensitivities
    from .prmg import mechanism_importance_values

    logger.info("Parsing RMG mechanism JSON file")

    with open(rmg_mech_json) as fle:
        mech_rxn_dcts = json.load(fle)

    sids = mechanism_species_identifiers(mech_rxn_dcts)
    rids = mechanism_reaction_identifiers(mech_rxn_dcts)
    ucrts = mechanism_uncertainties(mech_rxn_dcts)
    stvts = mechanism_sensitivities(mech_rxn_dcts)
    ipvls = mechanism_importance_values(mech_rxn_dcts)

    spc_df = pandas.DataFrame({'species_id': sids})
    rxn_df = pandas.DataFrame({'reaction_id': rids,
                               'uncertainty': ucrts,
                               'sensitivity': stvts,
                               'rmg_value': ipvls})

    logger.info("Writing species to {:s}".format(spc_csv_out))
    spc_df.to_csv(spc_csv_out, index=False)
    logger.info("Writing reactions to {:s}".format(rxn_csv_out))
    rxn_df.to_csv(rxn_csv_out, index=False)


def init_from_chemkin(chemkin_mech_txt, spc_csv, spc_csv_out, rxn_csv_out,
                      logger):
    """ initialize the mechanism from a CHEMKIN file
    """
    from .pchemkin import mechanism_reaction_identifiers

    logger.info("Parsing CHEMKIN mechanism text file")

    mech_str = open(chemkin_mech_txt).read()
    spc_df = pandas.read_csv(spc_csv)

    sid_dct = dict(zip(spc_df['species'], spc_df['species_id']))

    pi_rids, lp_rids, fo_rids = mechanism_reaction_identifiers(mech_str,
                                                               sid_dct)

    pi_df = pandas.DataFrame({'reaction_id': pi_rids})
    lp_df = pandas.DataFrame({'reaction_id': lp_rids})
    fo_df = pandas.DataFrame({'reaction_id': fo_rids})
    rxn_df = pandas.concat([pi_df, lp_df, fo_df],
                           keys=['press_indep', 'low_p', 'falloff'],
                           names=['pressure_dependence', ''])

    logger.info("Writing species to {:s}".format(spc_csv_out))
    spc_df.to_csv(spc_csv_out, index=False)
    logger.info("Writing reactions to {:s}".format(rxn_csv_out))
    rxn_df.to_csv(rxn_csv_out, index=False)


def abstractions_init(spc_csv, rxn_csv, rxn_csv_out, can_csv_out, logger,
                      keep_rxn_csv):
    """ initialize abstractions
    """
    logger.info("Reading in {:s}".format(rxn_csv))
    rxn_df = pandas.read_csv(rxn_csv)

    logger.info("Determining abstraction candidates by formula")
    rxn_df['abstraction_candidate'] = map(
        abstraction_candidate, rxn_df['reaction_id'])

    logger.info("Reading in species geometries from {:s}".format(spc_csv))
    mgeo_dct = geometries(spc_csv)

    logger.info("Iterating over candidates")
    abs_rows = []
    can_rows = []
    candidates = rxn_df['abstraction_candidate']
    for idx, rid in rxn_df[candidates]['reaction_id'].iteritems():
        logger.info('reaction {:d}: {:s}'.format(idx, rid))
        abstr = abstraction_indices(rid, mgeo_dct)
        if abstr:
            new_rid, idxs = abstr

            logger.info('  found abstraction!')
            rxn_df.loc[idx, 'class'] = 'abstraction'

            logger.info('  new reaction ID: {:s}'.format(new_rid))
            rxn_df.loc[idx, 'reaction_id'] = new_rid

            logger.info('  indices: {:s}'.format(str(idxs)))
            abs_rows.append((new_rid,) + idxs)
        else:
            can_rows.append((rid,))

    logger.info("Writing abstractions to {:s}"
                .format(os.path.abspath(rxn_csv_out)))
    rxn_df_out = pandas.DataFrame(
        abs_rows,
        columns=('reaction_id', 'q1h_idx', 'q2_idx', 'q1_idx', 'q2h_idx'))
    rxn_df_out.to_csv(rxn_csv_out, index=False)

    logger.info("Writing left-over candidates to {:s}"
                .format(os.path.abspath(can_csv_out)))
    can_df_out = pandas.DataFrame(
        can_rows,
        columns=('reaction_id',))
    can_df_out.to_csv(can_csv_out, index=False)

    rxn_csv = timestamp_if_exists(rxn_csv) if keep_rxn_csv else rxn_csv
    logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
    rxn_df.to_csv(rxn_csv, index=False)


def additions_init(spc_csv, rxn_csv, rxn_csv_out, can_csv_out, logger,
                   keep_rxn_csv):
    """ initialize additions
    """
    logger.info("Reading in {:s}".format(rxn_csv))
    rxn_df = pandas.read_csv(rxn_csv)

    logger.info("Determining addition candidates by formula")
    rxn_df['addition_candidate'] = map(
        addition_candidate, rxn_df['reaction_id'])

    logger.info("Reading in species geometries from {:s}".format(spc_csv))
    mgeo_dct = geometries(spc_csv)

    logger.info("Iterating over candidates")
    adt_rows = []
    can_rows = []
    candidates = rxn_df['addition_candidate']
    for idx, rid in rxn_df[candidates]['reaction_id'].iteritems():
        logger.info('reaction {:d}: {:s}'.format(idx, rid))
        addtn = addition_indices(rid, mgeo_dct)
        if addtn:
            new_rid, idxs = addtn

            logger.info('  found addition!')
            rxn_df.loc[idx, 'class'] = 'addition'

            logger.info('  new reaction ID: {:s}'.format(new_rid))
            rxn_df.loc[idx, 'reaction_id'] = new_rid

            logger.info('  indices: {:s}'.format(str(idxs)))
            adt_rows.append((new_rid,) + idxs)
        else:
            can_rows.append((rid,))

    logger.info("Writing additions to {:s}"
                .format(os.path.abspath(rxn_csv_out)))
    rxn_df_out = pandas.DataFrame(
        adt_rows,
        columns=('reaction_id', 'x_idx', 'y_idx', 'xy_idx_x', 'xy_idx_y'))
    rxn_df_out.to_csv(rxn_csv_out, index=False)

    logger.info("Writing left-over candidates to {:s}"
                .format(os.path.abspath(can_csv_out)))
    can_df_out = pandas.DataFrame(
        can_rows,
        columns=('reaction_id',))
    can_df_out.to_csv(can_csv_out, index=False)

    rxn_csv = timestamp_if_exists(rxn_csv) if keep_rxn_csv else rxn_csv
    logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
    rxn_df.to_csv(rxn_csv, index=False)


def abstractions_run(spc_csv, batch_csv, rxn_csv, tmp_txt, tmp_keyval_str, run_dir,
                     escape_id_chars, job_argv, logger, keep_rxn_csv):
    """ run abstractions
    """
    logger.info("Reading in {:s}".format(rxn_csv))
    rxn_df = pandas.read_csv(rxn_csv)

    logger.info("Reading in {:s}".format(batch_csv))
    batch_df = pandas.read_csv(batch_csv)

    logger.info("Reading in species geometries from {:s}".format(spc_csv))
    mgeo_dct = geometries(spc_csv)

    logger.info("Reading template file from {:s}".format(tmp_txt))
    tmp_str = read_file(tmp_txt)

    tmp_keyval_dct = _interpret_template_key_values(tmp_keyval_str)
    logger.info("Substitution key values: {:s}".format(str(tmp_keyval_dct)))

    if not os.path.exists(run_dir):
        logger.info("Creating run directory {:s}".format(run_dir))
        os.mkdir(run_dir)

    logger.info("Writing files for reactions in {:s}".format(batch_csv))
    batch_rids = tuple(batch_df['reaction_id'])
    select = list(map(batch_rids.__contains__, rxn_df['reaction_id']))
    for idx, row in rxn_df[select].iterrows():
        rid = row['reaction_id']
        logger.info('reaction {:d}: {:s}'.format(idx, rid))

        idxs = tuple(row[['q1h_idx', 'q2_idx']])
        logger.info('  q1h_idx={:d}, q2_idx={:d}'.format(*idxs))
        dxyz_dct = abstraction_xyz_string_dictionary(rid, idxs, mgeo_dct)
        logger.info(dxyz_dct)
        if dxyz_dct:
            sids, dxyzs = zip(*dxyz_dct.items())

            dname = escape_id_chars(rid)
            dpath = os.path.join(run_dir, dname)
            logger.info("Creating job directory {:s}".format(dpath))
            if not os.path.exists(dpath):
                os.mkdir(dpath)

            fnames = tuple(map('{:s}.xyz'.format, map(escape_id_chars, sids)))
            fpaths = tuple(os.path.join(dpath, fname) for fname in fnames)
            for fpath, dxyz in zip(fpaths, dxyzs):
                logger.info("  Writing {:s}".format(fpath))
                write_file(fpath, dxyz)

            inp_str = abstraction_input_string(rid, tmp_str, tmp_keyval_dct)
            inp_fpath = os.path.join(dpath, 'input.dat')

            logger.info("  Writing {:s}".format(inp_fpath))
            write_file(inp_fpath, inp_str)
            rxn_df.loc[idx, 'created'] = True
            rxn_df.loc[idx, 'path'] = dpath
        else:
            logger.info("  Failed to create .xyz files")
            rxn_df.loc[idx, 'created'] = False

    rxn_csv = timestamp_if_exists(rxn_csv) if keep_rxn_csv else rxn_csv
    logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
    rxn_df.to_csv(rxn_csv, index=False)

    rxn_csv = timestamp_if_exists(rxn_csv) if keep_rxn_csv else rxn_csv
    logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
    rxn_df.to_csv(rxn_csv, index=False)

    owd = os.getcwd()
    logger.info("Running job command in successfully created directories")
    for idx, row in rxn_df[select].iterrows():
        if row['created']:
            rid = row['reaction_id']
            logger.info('reaction {:d}: {:s}'.format(idx, rid))

            path = row['path']
            logger.info('  entering {:s}'.format(path))
            os.chdir(path)

            cmd_str = ' '.join(job_argv)
            logger.info("  running command '{:s}' in {:s}"
                        .format(cmd_str, path))
            try:
                subprocess.check_call(job_argv)
            except Exception as err:
                logger.info("   command '{:s}' failed with error '{:s}'"
                            .format(cmd_str, err))

            os.chdir(owd)


def additions_run(spc_csv, batch_csv, rxn_csv, tmp_txt, tmp_keyval_str, run_dir,
                  escape_id_chars, job_argv, logger, keep_rxn_csv):
    """ run additions
    """
    logger.info("Reading in {:s}".format(rxn_csv))
    rxn_df = pandas.read_csv(rxn_csv)

    logger.info("Reading in {:s}".format(batch_csv))
    batch_df = pandas.read_csv(batch_csv)

    logger.info("Reading in species geometries from {:s}".format(spc_csv))
    mgeo_dct = geometries(spc_csv)

    logger.info("Reading template file from {:s}".format(tmp_txt))
    tmp_str = read_file(tmp_txt)

    tmp_keyval_dct = _interpret_template_key_values(tmp_keyval_str)
    logger.info("Substitution key values: {:s}".format(str(tmp_keyval_dct)))

    if not os.path.exists(run_dir):
        logger.info("Creating run directory {:s}".format(run_dir))
        os.mkdir(run_dir)

    logger.info("Writing files for reactions in {:s}".format(batch_csv))
    batch_rids = tuple(batch_df['reaction_id'])
    select = list(map(batch_rids.__contains__, rxn_df['reaction_id']))
    for idx, row in rxn_df[select].iterrows():
        rid = row['reaction_id']
        logger.info('reaction {:d}: {:s}'.format(idx, rid))

        idxs = tuple(row[['x_idx', 'y_idx']])
        logger.info('  x_idx={:d}, y_idx={:d}'.format(*idxs))
        dxyz_dct = addition_xyz_string_dictionary(rid, idxs, mgeo_dct)
        if dxyz_dct:
            sids, dxyzs = zip(*dxyz_dct.items())

            dname = escape_id_chars(rid)
            dpath = os.path.join(run_dir, dname)
            logger.info("Creating job directory {:s}".format(dpath))
            if not os.path.exists(dpath):
                os.mkdir(dpath)

            fnames = tuple(map('{:s}.xyz'.format, map(escape_id_chars, sids)))
            fpaths = tuple(os.path.join(dpath, fname) for fname in fnames)
            for fpath, dxyz in zip(fpaths, dxyzs):
                logger.info("  Writing {:s}".format(fpath))
                write_file(fpath, dxyz)

            inp_str = addition_input_string(rid, tmp_str, tmp_keyval_dct)
            inp_fpath = os.path.join(dpath, 'input.dat')

            logger.info("  Writing {:s}".format(inp_fpath))
            write_file(inp_fpath, inp_str)
            rxn_df.loc[idx, 'created'] = True
            rxn_df.loc[idx, 'path'] = dpath
        else:
            logger.info("  Failed to create .xyz files")
            rxn_df.loc[idx, 'created'] = False

    rxn_csv = timestamp_if_exists(rxn_csv) if keep_rxn_csv else rxn_csv
    logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
    rxn_df.to_csv(rxn_csv, index=False)

    owd = os.getcwd()
    logger.info("Running job command in successfully created directories")
    for idx, row in rxn_df[select].iterrows():
        if row['created']:
            rid = row['reaction_id']
            logger.info('reaction {:d}: {:s}'.format(idx, rid))

            path = row['path']
            logger.info('  entering {:s}'.format(path))
            os.chdir(path)

            cmd_str = ' '.join(job_argv)
            logger.info("  running command '{:s}' in {:s}"
                        .format(cmd_str, path))
            try:
                subprocess.check_call(job_argv)
            except Exception as err:
                logger.info("   command '{:s}' failed with error '{:s}'"
                            .format(cmd_str, err))

            os.chdir(owd)


def _interpret_template_key_values(tmp_keyval_str):
    tmp_keyval_dct = dict(
        (s.strip() for s in kv.split(':')) for kv in tmp_keyval_str.split('|'))
    return tmp_keyval_dct
