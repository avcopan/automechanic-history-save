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
from .iohelp2 import abstraction_candidate
from .iohelp2 import abstraction
from .iohelp2 import abstraction_xyz_strings
from .iohelp2 import abstraction_input_string
from .iohelp2 import addition_candidate
from .iohelp2 import addition
from .iohelp2 import addition_xyz_strings
from .iohelp2 import addition_input_string
from .iohelp2 import migration_candidate
from .iohelp2 import migration
from .iohelp2 import migration_xyz_strings
from .iohelp2 import migration_input_string


# logging functions
def init(mech_txt, spc_csv, spc_csv_out, rxn_csv_out, geom_dir, id2path,
         logger):
    """ initialize a mechanism from a CHEMKIN mechanism file
    """
    from .iohelp2 import translate_chemkin_reaction
    from .iohelp2 import translate_chemkin_thermo_data
    from .pchemkin2 import reactions as chemkin_reactions
    from .pchemkin2 import therm_datas as chemkin_therm_datas

    logger.info("Reading in {:s}".format(mech_txt))
    mech_str = open(mech_txt).read()

    logger.info("Reading in {:s}".format(spc_csv))
    spc_df = pandas.read_csv(spc_csv)

    spc_df['species_id'] = map(canonical_species_identifier,
                               spc_df['species_id'])

    sid_dct = dict(zip(spc_df['species'], spc_df['species_id']))

    logger.info("Finding reactions")
    rxn_strs = chemkin_reactions(mech_str)
    rxn_rows = []
    mis_rows = []
    for num, rxn_str in enumerate(rxn_strs):
        rid = translate_chemkin_reaction(rxn_str, sid_dct)
        if rid:
            rid = canonical_reaction_identifier(rid)
            logger.info("Found reaction {:s}".format(rid))
            rxn_rows.append((rid, rxn_str, num))
        else:
            logger.info("Failed to translate reaction {:s}".format(rxn_str))
            mis_rows.append((rxn_str, num))

    thd_strs = chemkin_therm_datas(mech_str)
    thv_dct = dict(filter(bool,
                          (translate_chemkin_thermo_data(thd_str, sid_dct)
                           for thd_str in thd_strs)))

    spc_df = initialize_thermo_data(spc_df, thv_dct, logger)
    spc_df = initialize_geometries(spc_df, geom_dir, id2path, logger)

    rxn_cols = ('reaction_id', 'reaction', 'number')
    rxn_df = pandas.DataFrame(rxn_rows, columns=rxn_cols)

    mis_cols = ('reaction', 'number')
    mis_df = pandas.DataFrame(mis_rows, columns=mis_cols)

    logger.info("Writing species to {:s}".format(spc_csv_out))
    spc_df.to_csv(spc_csv_out, index=False)

    logger.info("Writing reactions to {:s}".format(rxn_csv_out))
    rxn_df.to_csv(rxn_csv_out, index=False)

    logger.info("Writing missed reactions to {:s}".format(rxn_csv_out))
    mis_df.to_csv('missed.csv', index=False)


def init_from_rmg(rmg_mech_json, spc_csv_out, rxn_csv_out, geom_dir, id2path,
                  logger):
    """ initialize a mechanism from RMG's JSON file
    """
    from .prmg import mechanism_species_identifiers
    from .prmg import mechanism_reaction_identifiers
    from .prmg import mechanism_sensitivities
    from .prmg import mechanism_importance_values

    logger.info("Parsing RMG mechanism JSON file")

    with open(rmg_mech_json) as fle:
        mech_rxn_dcts = json.load(fle)

    sids = mechanism_species_identifiers(mech_rxn_dcts)
    rids = mechanism_reaction_identifiers(mech_rxn_dcts)
    stvts = mechanism_sensitivities(mech_rxn_dcts)
    ipvls = mechanism_importance_values(mech_rxn_dcts)

    spc_df = pandas.DataFrame({'species_id': sids})
    rxn_df = pandas.DataFrame({'reaction_id': rids,
                               'sensitivity': stvts,
                               'rmg_value': ipvls})

    logger.info("Canonicalizing species IDs")
    spc_df['species_id'] = map(canonical_species_identifier,
                               spc_df['species_id'])

    logger.info("Canonicalizing reaction IDs")
    rxn_df['reaction_id'] = map(canonical_reaction_identifier,
                                rxn_df['reaction_id'])

    spc_df = initialize_geometries(spc_df, geom_dir, id2path, logger)

    logger.info("Writing species to {:s}".format(spc_csv_out))
    spc_df.to_csv(spc_csv_out, index=False)

    logger.info("Writing reactions to {:s}".format(rxn_csv_out))
    rxn_df.to_csv(rxn_csv_out, index=False)


def initialize_thermo_data(spc_df, thv_dct, logger):
    """ add thermo data to the species table
    """
    assert 'species_id' in spc_df

    for idx, row in spc_df.iterrows():
        sid = row['species_id']
        if sid in thv_dct:
            thv = thv_dct[sid]
            logger.info("Found thermo data for {:s}, H298 = {:f}"
                        .format(sid, thv))
            spc_df.at[idx, 'therm_val'] = thv
        else:
            logger.info("No thermo data for {:s}".format(sid))

    return spc_df


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
    rxn_df1.to_csv(rxn_csv1_out)

    rxn_csv2_out = os.path.join(dir2, rxn_csv_out)
    logger.info("Writing out-of-category reactions to {:s}"
                .format(rxn_csv2_out))
    if not os.path.exists(dir2):
        os.mkdir(dir2)
    rxn_df2.to_csv(rxn_csv2_out)

    logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
    timestamp_if_exists(rxn_csv)
    rxn_df.to_csv(rxn_csv, index=False)


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
        cols = (('reaction_id',) + sid_cols + idx_cols)
        rxn_df_out = pandas.DataFrame(rxn_rows, columns=cols)
        rxn_df_out.to_csv(rxn_csv_out, index=False)

        logger.info("Writing left-over candidates to {:s}"
                    .format(os.path.abspath(cdt_csv_out)))
        columns = ('reaction_id', 'exception')
        cdt_df_out = pandas.DataFrame(cdt_rows, columns=columns)
        cdt_df_out.to_csv(cdt_csv_out, index=False)

        logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
        timestamp_if_exists(rxn_csv)
        rxn_df.to_csv(rxn_csv, index=False)

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

            logger.info('here -1')
            dxyz_dct = reaction_xyz_strings(sids, idxs, mgeo_dct)
            if dxyz_dct:
                logger.info('here 0')
                dxyz_sids = dxyz_dct.keys()
                dxyzs = dxyz_dct.values()
                logger.info('here 1')

                dname = id2path(rid)
                logger.info('here 2')
                dpath = os.path.join(run_dir, dname)
                logger.info("Creating job directory {:s}".format(dpath))
                if not os.path.exists(dpath):
                    os.mkdir(dpath)

                logger.info('here 3')
                fnames = tuple(map('{:s}.xyz'.format, map(id2path, dxyz_sids)))
                logger.info('here 4')
                fpaths = tuple(os.path.join(dpath, fname) for fname in fnames)
                for fpath, dxyz in zip(fpaths, dxyzs):
                    logger.info("  Writing {:s}".format(fpath))
                    write_file(fpath, dxyz)
                logger.info('here 5')

                inp_str = reaction_input_string(sids, tmp_str, tmp_keyval_dct)
                inp_fpath = os.path.join(dpath, 'input.dat')
                logger.info('here 6')

                logger.info('here 7')
                logger.info("  Writing {:s}".format(inp_fpath))
                write_file(inp_fpath, inp_str)
                rxn_df.loc[idx, 'created'] = True
                rxn_df.loc[idx, 'path'] = dpath
            else:
                logger.info("  Failed to create .xyz files")
                rxn_df.loc[idx, 'created'] = False

        logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
        timestamp_if_exists(rxn_csv)
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
