""" drivers and functions for I/O
"""
import os
import time
import json
import subprocess
from functools import partial
import pandas
from .strid import canonical as canonical_species_identifier
from .strid import canonical_reaction_identifier
from .iohelp import abstraction_candidate
from .iohelp import abstraction
from .iohelp import abstraction_xyz_strings
from .iohelp import abstraction_input_string
from .iohelp import addition_candidate
from .iohelp import addition
from .iohelp import addition_xyz_strings
from .iohelp import addition_input_string
from .iohelp import migration_candidate
from .iohelp import migration
from .iohelp import migration_xyz_strings
from .iohelp import migration_input_string
from .table import from_columns as table_from_columns
from .table import from_rows as table_from_rows
from .table import reindex as reindex_table
from .table import sort as sort_table
from .table import merge as merge_tables


def init(mech_txt, spc_csv, rxn_csv_out, spc_csv_out, geom_dir, id2path,
         therm_txt, without_thermo, logger):
    """ initialize a mechanism from a CHEMKIN mechanism file
    """
    from .iohelp import translate_chemkin_reaction
    from .iohelp import thermo_value_dictionary
    from .pchemkin import reactions as chemkin_reactions
    from .pchemkin import therm_data_strings as chemkin_therm_data_strings

    logger.info("Reading in {:s}".format(mech_txt))
    mech_str = read_file(mech_txt)

    logger.info("Reading in {:s}".format(spc_csv))
    spc_df = pandas.read_csv(spc_csv)

    spcs = tuple(spc_df['species'])
    sids = tuple(map(canonical_species_identifier, spc_df['species_id']))
    sid_dct = dict(zip(spcs, sids))

    spc_df['species_id'] = sids

    logger.info("Finding reactions")
    rxn_strs = chemkin_reactions(mech_str)
    rxn_rows = []
    mis_rows = []
    for num, rxn_str in enumerate(rxn_strs):
        rid = translate_chemkin_reaction(rxn_str, sid_dct)
        if rid:
            rid = canonical_reaction_identifier(rid)
            logger.info("Found reaction {:s}".format(rid))
            rxn_rows.append((rid, num+1, rxn_str))
        else:
            logger.info("Failed to translate reaction {:s}".format(rxn_str))
            mis_rows.append((num+1, rxn_str))

    if not without_thermo:
        if therm_txt is None:
            therm_str = mech_str
        else:
            logger.info("Reading in {:s}".format(therm_txt))
            therm_str = read_file(therm_txt)
        thd_strs = chemkin_therm_data_strings(therm_str)
        thv_dct = thermo_value_dictionary(thd_strs, sid_dct)
        spc_df['therm_val'] = map(thv_dct.__getitem__, spc_df['species_id'])

    spc_df = initialize_geometries(spc_df, geom_dir, id2path, logger)

    rxn_cols = ('reaction_id', 'chemkin_index', 'reaction')
    rxn_df = table_from_rows(rxn_rows, rxn_cols)

    mis_cols = ('chemkin_index', 'reaction')
    mis_df = table_from_rows(mis_rows, mis_cols)

    logger.info("Writing species to {:s}".format(spc_csv_out))
    write_table_to_csv(spc_df, spc_csv_out, float_fmt='%.8f')

    logger.info("Writing reactions to {:s}".format(rxn_csv_out))
    write_table_to_csv(rxn_df, rxn_csv_out, float_fmt='%.4f')

    logger.info("Writing missed reactions to {:s}".format(rxn_csv_out))
    write_table_to_csv(mis_df, 'missed.csv')


def init_from_rmg(mech_json, spc_json, rxn_csv_out, spc_csv_out, geom_dir,
                  id2path, logger):
    """ initialize a mechanism from RMG files
    """
    from .prmg import species_name as species_name_from_dct
    from .prmg import species_identifier as species_identifier_from_dct
    from .prmg import species_thermo_value as species_thermo_value_from_dct
    from .prmg import reaction_name as reaction_name_from_dct
    from .prmg import reaction_identifier as reaction_identifier_from_dct
    from .prmg import reaction_sensitivity as reaction_sensitivity_from_dct
    from .prmg import reaction_uncertainty as reaction_uncertainty_from_dct
    from .prmg import reaction_value as reaction_value_from_dct

    logger.info("Reading in {:s}".format(mech_json))
    mech_rxn_dcts = read_json(mech_json)

    logger.info("Reading in {:s}".format(spc_json))
    spc_dcts = read_json(spc_json)

    spc_strs = list(map(species_name_from_dct, spc_dcts))
    spc_sids = list(map(canonical_species_identifier,
                        map(species_identifier_from_dct, spc_dcts)))
    spc_thvs = list(map(species_thermo_value_from_dct, spc_dcts))

    mech_rxn_strs = list(map(reaction_name_from_dct, mech_rxn_dcts))
    mech_rids = list(map(canonical_reaction_identifier,
                         map(reaction_identifier_from_dct, mech_rxn_dcts)))
    mech_stvys = list(map(reaction_sensitivity_from_dct, mech_rxn_dcts))
    mech_uctys = list(map(reaction_uncertainty_from_dct, mech_rxn_dcts))
    mech_rvals = list(map(reaction_value_from_dct, mech_rxn_dcts))

    spc_cols = (spc_sids, spc_strs, spc_thvs)
    spc_col_keys = ('species_id', 'species', 'therm_val')
    spc_df = table_from_columns(spc_cols, spc_col_keys)
    spc_df = initialize_geometries(spc_df, geom_dir, id2path, logger)

    rxn_cols = (mech_rids, mech_rxn_strs, mech_stvys, mech_uctys, mech_rvals)
    rxn_col_keys = ('reaction_id', 'reaction', 'sensitivity', 'uncertainty',
                    'rmg_value')
    rxn_df = table_from_columns(rxn_cols, rxn_col_keys)

    logger.info("Writing species to {:s}".format(spc_csv_out))
    write_table_to_csv(spc_df, spc_csv_out, float_fmt='%.8f')

    logger.info("Writing reactions to {:s}".format(rxn_csv_out))
    write_table_to_csv(rxn_df, rxn_csv_out, float_fmt='%.4f')


def initialize_geometries(spc_df, geom_dir, id2path, logger):
    """ initialize species geometries
    """
    from .strid import xyz_string

    logger.info("Initializing species geometries")

    assert 'species_id' in spc_df

    if not os.path.exists(geom_dir):
        os.mkdir(geom_dir)

    for idx, row in spc_df.iterrows():
        sid = row['species_id']
        fname = id2path(sid) + '.xyz'
        fpath = os.path.join(geom_dir, fname)
        if not os.path.exists(fpath):
            logger.info("Writing geometry for {:s} to {:s}"
                        .format(sid, fpath))
            dxyz = xyz_string(sid)
            write_file(fpath, contents=dxyz)
        else:
            logger.info("Geometry for {:s} already exists at {:s}"
                        .format(sid, fpath))
        spc_df.at[idx, 'path'] = fpath

    return spc_df


def read_geometries(spc_csv):
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


def read_thermo_data(spc_csv):
    """ a dictionary of thermo values (H298), indexed by species ID
    """
    thv_dct = None

    spc_df = pandas.read_csv(spc_csv)
    if 'therm_val' in spc_df:
        thv_dct = dict(zip(spc_df['species_id'], spc_df['therm_val']))

    return thv_dct


def abstractions_init(spc_csv, rxn_csv, rxn_csv_out, cdt_csv_out, logger):
    """ initialize abstraction reactions
    """
    _init = reactions_initializer(
        cls='abstraction',
        is_candidate=abstraction_candidate,
        reaction=abstraction,
        sid_cols=('q1h', 'q2', 'q1', 'q2h'),
        idx_cols=('q1h_idx', 'q2_idx', 'q1_idx', 'q2h_idx')
    )
    return _init(spc_csv, rxn_csv, rxn_csv_out, cdt_csv_out, logger)


def additions_init(spc_csv, rxn_csv, rxn_csv_out, cdt_csv_out, logger):
    """ initialize addition reactions
    """
    _init = reactions_initializer(
        cls='addition',
        is_candidate=addition_candidate,
        reaction=addition,
        sid_cols=('x', 'y', 'xy'),
        idx_cols=('x_idx', 'y_idx', 'xy_idx_x', 'xy_idx_y')
    )
    return _init(spc_csv, rxn_csv, rxn_csv_out, cdt_csv_out, logger)


def migrations_init(spc_csv, rxn_csv, rxn_csv_out, cdt_csv_out, logger):
    """ initialize migration reactions
    """
    _init = reactions_initializer(
        cls='migration',
        is_candidate=migration_candidate,
        reaction=migration,
        sid_cols=('r', 'p'),
        idx_cols=('r_idx_h', 'r_idx_a', 'p_idx_h', 'p_idx_a')
    )
    return _init(spc_csv, rxn_csv, rxn_csv_out, cdt_csv_out, logger)


def abstractions_run(spc_csv, batch_csv, rxn_csv, tmp_txt, tmp_keyval_str,
                     run_dir, id2path, job_argv, logger):
    """ run abstractions
    """
    _run = reactions_runner(
        cls='abstraction',
        reaction_xyz_strings=abstraction_xyz_strings,
        reaction_input_string=abstraction_input_string,
        sid_cols=('q1h', 'q2', 'q1', 'q2h'),
        idx_cols=('q1h_idx', 'q2_idx', 'q1_idx', 'q2h_idx')
    )
    return _run(spc_csv, batch_csv, rxn_csv, tmp_txt, tmp_keyval_str, run_dir,
                id2path, job_argv, logger)


def additions_run(spc_csv, batch_csv, rxn_csv, tmp_txt, tmp_keyval_str,
                  run_dir, id2path, job_argv, logger):
    """ run additions
    """
    _run = reactions_runner(
        cls='addition',
        reaction_xyz_strings=addition_xyz_strings,
        reaction_input_string=addition_input_string,
        sid_cols=('x', 'y', 'xy'),
        idx_cols=('x_idx', 'y_idx', 'xy_idx_x', 'xy_idx_y')
    )
    return _run(spc_csv, batch_csv, rxn_csv, tmp_txt, tmp_keyval_str, run_dir,
                id2path, job_argv, logger)


def migrations_run(spc_csv, batch_csv, rxn_csv, tmp_txt, tmp_keyval_str,
                   run_dir, id2path, job_argv, logger):
    """ run migrations
    """
    _run = reactions_runner(
        cls='migration',
        reaction_xyz_strings=migration_xyz_strings,
        reaction_input_string=migration_input_string,
        sid_cols=('r', 'p'),
        idx_cols=('r_idx_h', 'r_idx_a', 'p_idx_h', 'p_idx_a')
    )
    return _run(spc_csv, batch_csv, rxn_csv, tmp_txt, tmp_keyval_str, run_dir,
                id2path, job_argv, logger)


def divide(key, dir1, dir2, rxn_csv, rxn_csv_out, logger):
    """ split reactions by key
    """
    from .strid import is_radical_radical
    from .strid import is_spin_balanced

    if key == 'rad-rad':
        meets_condition_ = is_radical_radical
    elif key == 'high-spin':
        meets_condition_ = is_spin_balanced
    else:
        raise ValueError("Unrecognized divide key: {:s}".format(key))

    logger.info("Reading in {:s}".format(rxn_csv))
    rxn_df = pandas.read_csv(rxn_csv)
    rxn_df[key] = map(meets_condition_, rxn_df['reaction_id'])
    rxn_df1 = rxn_df[rxn_df[key]].drop(columns=key)
    rxn_df2 = rxn_df[~rxn_df[key]].drop(columns=key)

    rxn_csv1_out = os.path.join(dir1, rxn_csv_out)
    logger.info("Writing in-category reactions to {:s}"
                .format(rxn_csv1_out))
    if not os.path.exists(dir1):
        os.mkdir(dir1)
    write_table_to_csv(rxn_df1, rxn_csv1_out)

    rxn_csv2_out = os.path.join(dir2, rxn_csv_out)
    logger.info("Writing out-of-category reactions to {:s}"
                .format(rxn_csv2_out))
    if not os.path.exists(dir2):
        os.mkdir(dir2)
    write_table_to_csv(rxn_df2, rxn_csv2_out)

    logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
    write_table_to_csv(rxn_df, rxn_csv)


def csv_reindex(table_csv, logger):
    """ reindex a table
    """
    logger.info("Reading in {:s}".format(table_csv))
    table_df = pandas.read_csv(table_csv)

    table_df = reindex_table(table_df)

    logger.info("Writing updated {:s}".format(table_csv))
    write_table_to_csv(table_df, table_csv)


def csv_sort(table_csv, col_key, descending, logger):
    """ sort table by column
    """
    logger.info("Reading in {:s}".format(table_csv))
    table_df = pandas.read_csv(table_csv)

    table_df = sort_table(table_df, col_key, descending=descending)

    logger.info("Writing updated {:s}".format(table_csv))
    write_table_to_csv(table_df, table_csv)


def csv_merge(table_csvs, col_key, table_csv_out, logger):
    """ merge tables by column
    """
    table_dfs = []
    for table_csv in table_csvs:
        logger.info("Reading in {:s}".format(table_csv))
        table_df = pandas.read_csv(table_csv)
        table_dfs.append(table_df)

    table_df_out = merge_tables(table_dfs, col_key)

    logger.info("Writing {:s}".format(table_csv_out))
    write_table_to_csv(table_df_out, table_csv_out)


# meta scripts
def reactions_initializer(cls, is_candidate, reaction, sid_cols, idx_cols):
    """ initialize reactions
    """
    assert cls in ('abstraction', 'addition', 'migration')

    def _init(spc_csv, rxn_csv, rxn_csv_out, cdt_csv_out, logger):
        logger.info("Reading in {:s}".format(rxn_csv))
        rxn_df = pandas.read_csv(rxn_csv)

        logger.info("Reading in species geometries from {:s}".format(spc_csv))
        mgeo_dct = read_geometries(spc_csv)

        logger.info("Reading thermo data from {:s}".format(spc_csv))
        thv_dct = read_thermo_data(spc_csv)
        if not thv_dct:
            logger.info("No thermo data found.")

        logger.info("Iterating over candidates")
        rxn_rows = []
        cdt_rows = []
        for idx, rid in rxn_df['reaction_id'].iteritems():
            if is_candidate(rid):
                logger.info('reaction {:d}: {:s}'.format(idx, rid))

                err = None
                try:
                    rxn = reaction(rid, mgeo_dct, thv_dct)
                except Exception as err:
                    logger.info('  exception: {:s}!'.format(str(err)))
                    rxn = None

                if rxn:
                    sids, idxs = rxn

                    logger.info('  found {:s}!'.format(cls))
                    rxn_df.loc[idx, 'class'] = cls

                    log_sids = ', '.join(
                        '{:s}: {:s}'.format(sid_col, sid)
                        for sid_col, sid in zip(sid_cols, sids))
                    log_idxs = ', '.join(
                        '{:s}: {:d}'.format(idx_col, idx)
                        for idx_col, idx in zip(idx_cols, idxs))

                    logger.info('  {:s}\n  {:s}'.format(log_sids, log_idxs))

                    rxn_rows.append((rid,) + sids + idxs)
                else:
                    cdt_rows.append((rid, err))

        logger.info("Writing {:s} reactions to {:s}"
                    .format(cls, os.path.abspath(rxn_csv_out)))
        col_keys = (('reaction_id',) + sid_cols + idx_cols)
        rxn_df_out = table_from_rows(rxn_rows, col_keys)
        write_table_to_csv(rxn_df_out, rxn_csv_out)

        logger.info("Writing left-over candidates to {:s}"
                    .format(os.path.abspath(cdt_csv_out)))
        cdt_col_keys = ('reaction_id', 'exception')
        cdt_df_out = table_from_rows(cdt_rows, cdt_col_keys)
        write_table_to_csv(cdt_df_out, cdt_csv_out)

        logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
        write_table_to_csv(rxn_df, rxn_csv)

    return _init


def reactions_runner(cls, reaction_xyz_strings, reaction_input_string,
                     sid_cols, idx_cols):
    """ run reactions
    """
    assert cls in ('abstraction', 'addition', 'migration')

    def _run(spc_csv, batch_csv, rxn_csv, tmp_txt, tmp_keyval_str, run_dir,
             id2path, job_argv, logger):
        logger.info("Reading in {:s}".format(rxn_csv))
        rxn_df = pandas.read_csv(rxn_csv)

        logger.info("Reading in {:s}".format(batch_csv))
        batch_df = pandas.read_csv(batch_csv)

        logger.info("Reading in species geometries from {:s}".format(spc_csv))
        mgeo_dct = read_geometries(spc_csv)

        logger.info("Reading template file from {:s}".format(tmp_txt))
        tmp_str = read_file(tmp_txt)

        tmp_keyval_dct = _interpret_template_key_values(tmp_keyval_str)
        logger.info("Substitution key values: {:s}"
                    .format(str(tmp_keyval_dct)))

        if not os.path.exists(run_dir):
            logger.info("Creating run directory {:s}".format(run_dir))
            os.mkdir(run_dir)

        logger.info("Writing files for reactions in {:s}".format(batch_csv))
        batch_rids = tuple(batch_df['reaction_id'])
        select = list(map(batch_rids.__contains__, rxn_df['reaction_id']))
        for idx, row in rxn_df[select].iterrows():
            rid = row['reaction_id']
            logger.info('reaction {:d}: {:s}'.format(idx, rid))

            sids = tuple(row[list(sid_cols)])
            idxs = tuple(row[list(idx_cols)])
            logger.info('  indices: {:s}'.format(str(idxs)))

            dxyz_dct = reaction_xyz_strings(sids, idxs, mgeo_dct)
            if dxyz_dct:
                dxyz_sids = dxyz_dct.keys()
                dxyzs = dxyz_dct.values()

                dname = id2path(rid)
                dpath = os.path.join(run_dir, dname)
                logger.info("Creating job directory {:s}".format(dpath))
                if not os.path.exists(dpath):
                    os.mkdir(dpath)

                fnames = tuple(map('{:s}.xyz'.format, map(id2path, dxyz_sids)))
                fpaths = tuple(os.path.join(dpath, fname) for fname in fnames)
                for fpath, dxyz in zip(fpaths, dxyzs):
                    logger.info("  Writing {:s}".format(fpath))
                    write_file(fpath, dxyz)

                inp_str = reaction_input_string(sids, tmp_str, tmp_keyval_dct)
                inp_fpath = os.path.join(dpath, 'input.dat')

                logger.info("  Writing {:s}".format(inp_fpath))
                write_file(inp_fpath, inp_str)
                rxn_df.loc[idx, 'created'] = True
                rxn_df.loc[idx, 'path'] = dpath
            else:
                logger.info("  Failed to create .xyz files")
                rxn_df.loc[idx, 'created'] = False

        logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
        write_table_to_csv(rxn_df, rxn_csv)

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

    return _run


# non-logging functions
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


def read_json(fpath):
    """ read json file
    """
    return json.load(open(fpath))


def write_table_to_csv(table_df, table_csv, float_fmt=None):
    """ write table to csv
    """
    timestamp_if_exists(table_csv)
    table_df.to_csv(table_csv, index=False, float_format=float_fmt)


def timestamp_if_exists(fpath):
    """ open a file, avoiding overwrites if requested
    """
    if os.path.isfile(fpath):
        time_stamp = time.strftime("%Y%m%d-%H%M%S")
        new_fpath = "{:s}_{:s}".format(fpath, time_stamp)
        os.rename(fpath, new_fpath)


# helpers
def _interpret_template_key_values(tmp_keyval_str):
    tmp_keyval_dct = dict(
        (s.strip() for s in kv.split(':')) for kv in tmp_keyval_str.split('|'))
    return tmp_keyval_dct
