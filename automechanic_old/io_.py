""" drivers and functions for I/O
"""
from __future__ import unicode_literals
from builtins import open, str
import os
import time
import json
import subprocess
from itertools import chain
from itertools import starmap
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
from .table import set_column as set_table_column
from .table import is_empty_value as is_empty_table_value
from .table import has_column_keys as table_has_column_keys
from .table import column as table_column
from .table import columns as table_columns
from .table import lookup_dictionary as table_lookup_dictionary
from .table import lookup_row as table_lookup_row
from .table import lookup_update as table_lookup_update
from .table import columns_like as table_with_columns_like
from .table import update_column_keys as update_table_column_keys
from .table import append_rows as append_table_row_dicts
from .table import append_columns as append_table_columns
from .table import column_keys as table_column_keys
from .table import iterate_rows as iterate_table_row_dicts
from .table import from_columns as table_from_columns
from .table import from_rows as table_from_row_dicts
from .table import reindex as reindex_table
from .table import sort as sort_table
from .table import merge as merge_tables
from .table import intersect as intersect_tables
from .table import move_column_to_front as move_table_column_to_front
# from new table interface
from .table2 import EMPTY
from .table2 import row_indices
from .table2 import update_column_by_index
from .table2 import sql_where_eq
from .table2 import sql_where_in
from .table2 import sql_select_one
# from new parsing module
from .parse.rere.find import split
from .parse.rere.find import strip_spaces
# threading
from .thread import tag_team_starmap

ADD_XYZ_EXTENSION = '{:s}.xyz'.format
SID_COL_KEY = 'species_id'
GEOM_PATH_COL_KEY = 'geom_path'
RXN_IDX_COL_KEY = 'index'
RID_COL_KEY = 'reaction_id'
RXN_PATH_COL_KEY = 'path'
RXN_STAT_COL_KEY = 'status'
RXN_CREATED_VAL = 'created'
RXN_NOT_CREATED_VAL = 'not created'
RXN_RAN_VAL = 'ran'
RXN_FAILED_VAL = 'failed'
ARRH_COL_KEYS = ('arrh_a', 'arrh_b', 'arrh_e')
NASA_LO_COL_KEYS = ('nasa_lo_1', 'nasa_lo_2', 'nasa_lo_3', 'nasa_lo_4',
                    'nasa_lo_5', 'nasa_lo_6', 'nasa_lo_7')
NASA_HI_COL_KEYS = ('nasa_hi_1', 'nasa_hi_2', 'nasa_hi_3', 'nasa_hi_4',
                    'nasa_hi_5', 'nasa_hi_6', 'nasa_hi_7')
NASA_T_COL_KEYS = ('nasa_t_com', 'nasa_t_lo', 'nasa_t_hi')
REACTANT_SID_COL_KEYS = (
    ('addition', ('x', 'y')),
    ('abstraction', ('q1h', 'q2')),
    ('migration', ('r',))
)
PRODUCT_SID_COL_KEYS = (
    ('addition', ('xy')),
    ('abstraction', ('q1', 'q2h')),
    ('migration', ('p',))
)
REACTION_SID_COL_KEYS = (
    ('addition', ('x', 'y', 'xy')),
    ('abstraction', ('q1h', 'q2', 'q1', 'q2h')),
    ('migration', ('r', 'p'))
)
REACTION_IDX_COL_KEYS = (
    ('addition', ('x_idx', 'y_idx', 'xy_idx_x', 'xy_idx_y')),
    ('abstraction', ('q1h_idx', 'q2_idx', 'q1_idx', 'q2h_idx')),
    ('migration', ('r_idx_h', 'r_idx_a', 'p_idx_h', 'p_idx_a'))
)
REACTION_CANDIDATE_FINDERS = (
    ('addition', addition_candidate),
    ('abstraction', abstraction_candidate),
    ('migration', migration_candidate)
)
REACTION_FINDERS = (
    ('addition', addition),
    ('abstraction', abstraction),
    ('migration', migration)
)
REACTION_XYZ_STRING_MAKERS = (
    ('addition', addition_xyz_strings),
    ('abstraction', abstraction_xyz_strings),
    ('migration', migration_xyz_strings)
)
REACTION_INPUT_STRING_MAKERS = (
    ('addition', addition_input_string),
    ('abstraction', abstraction_input_string),
    ('migration', migration_input_string)
)


def init(mech_txt, spc_csv, rxn_csv_out, spc_csv_out, geom_dir, id2path,
         therm_txt, without_thermo, logger):
    """ initialize a mechanism from a CHEMKIN mechanism file
    """
    from .iohelp import translate_chemkin_reaction
    from .iohelp import thermo_value_dictionary
    from .pchemkin import reactions as chemkin_reactions
    from .pchemkin import thermo_block as chemkin_thermo_block
    from .pchemkin import therm_data_strings as chemkin_therm_data_strings

    logger.info("Reading in {:s}".format(mech_txt))
    mech_str = read_file(mech_txt)

    logger.info("Reading in {:s}".format(spc_csv))
    spc_df = pandas.read_csv(spc_csv)

    spcs = tuple(spc_df['species'])
    sids = tuple(map(canonical_species_identifier, spc_df['species_id']))
    sid_dct = dict(zip(spcs, sids))

    spc_df['species_id'] = sids

    if not without_thermo:
        if therm_txt is None and chemkin_thermo_block(mech_str):
            therm_str = mech_str
        else:
            logger.info("Reading in {:s}".format(therm_txt))
            therm_str = read_file(therm_txt)
        if not therm_str:
            raise ValueError("No thermo data found! Either specify the thermo "
                             "file or turn thermo off.")
        thd_strs = chemkin_therm_data_strings(therm_str)
        thv_dct = thermo_value_dictionary(thd_strs, sid_dct)
        spc_df['therm_val'] = tuple(
            map(thv_dct.__getitem__, spc_df['species_id']))

    logger.info("Finding reactions")
    rxn_strs = chemkin_reactions(mech_str)
    rxn_rows = []
    mis_rows = []
    seen = []
    for num, rxn_str in enumerate(rxn_strs):

        if rxn_str in seen:
            logger.info("Ignoring duplicate reaction {:s}".format(rxn_str))
            continue
        else:
            seen.append(rxn_str)

        ck_num = num + 1
        rid = translate_chemkin_reaction(rxn_str, sid_dct)
        if rid:
            rid = canonical_reaction_identifier(rid)
            logger.info("Found reaction {:d}: {:s}".format(ck_num, rid))
            rxn_rows.append((rid, ck_num, rxn_str))
        else:
            logger.info("Failed to translate reaction {:d}: {:s}"
                        .format(ck_num, rxn_str))
            mis_rows.append((ck_num, rxn_str))

    spc_df = initialize_geometries(spc_df, geom_dir, id2path, logger)

    rxn_cols = ('reaction_id', 'chemkin_index', 'reaction')
    rxn_df = table_from_row_dicts(rxn_rows, rxn_cols)
    rxn_df = sort_table(rxn_df, 'chemkin_index')

    mis_cols = ('chemkin_index', 'reaction')
    mis_df = table_from_row_dicts(mis_rows, mis_cols)

    logger.info("Writing species to {:s}".format(spc_csv_out))
    write_table_to_csv(spc_df, spc_csv_out, float_fmt='%.8f')

    logger.info("Writing reactions to {:s}".format(rxn_csv_out))
    write_table_to_csv(rxn_df, rxn_csv_out, float_fmt='%.4f')

    logger.info("Writing missed reactions to missed.csv")
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


def species_find_geometries(spc_csv, spc_csv_out, geom_dir, id2path, logger):
    """ find species .xyz files
    """

    logger.info("Reading in {:s}".format(spc_csv))
    spc_df = pandas.read_csv(spc_csv)

    spc_df = update_table_column_keys(spc_df, (GEOM_PATH_COL_KEY,))

    sids = table_column(spc_df, SID_COL_KEY)
    to_path_ = partial(os.path.join, geom_dir)
    fpaths = tuple(map(to_path_, map(ADD_XYZ_EXTENSION, map(id2path, sids))))

    for sid, fpath in zip(sids, fpaths):
        logger.info("species {:s}".format(sid))
        if os.path.exists(fpath):
            logger.info("  geometry file found at {:s}".format(fpath))
            spc_df = table_lookup_update(spc_df,
                                         (SID_COL_KEY, sid),
                                         (GEOM_PATH_COL_KEY, fpath))
        else:
            logger.info("  no geometry file found.")

    logger.info("Writing species to {:s}".format(spc_csv_out))
    write_table_to_csv(spc_df, spc_csv_out)


def species_fill_geometries(spc_csv, spc_csv_out, geom_dir, id2path, logger):
    """ find species .xyz files
    """
    from .strid import xyz_string

    logger.info("Reading in {:s}".format(spc_csv))
    spc_df = pandas.read_csv(spc_csv)

    if not os.path.exists(geom_dir):
        os.mkdir(geom_dir)

    spc_df = update_table_column_keys(spc_df, (GEOM_PATH_COL_KEY,))
    sids = table_column(spc_df, SID_COL_KEY)
    geom_fpaths = table_column(spc_df, GEOM_PATH_COL_KEY)

    pathless_sids = tuple(sid for sid, fpath in zip(sids, geom_fpaths)
                          if is_empty_table_value(fpath))

    to_fpath_ = partial(os.path.join, geom_dir)
    new_geom_fpaths = tuple(
        map(to_fpath_, map(ADD_XYZ_EXTENSION, map(id2path, pathless_sids))))

    for sid, fpath in zip(pathless_sids, new_geom_fpaths):
        logger.info("species {:s}".format(sid))
        if os.path.exists(fpath):
            logger.info("  file already exists; use geometry finder to add "
                        "  it to the database")
        else:
            logger.info("Writing geometry for {:s} to {:s}"
                        .format(sid, fpath))
            dxyz = xyz_string(sid)
            write_file(fpath, contents=dxyz)
            spc_df = table_lookup_update(spc_df,
                                         (SID_COL_KEY, sid),
                                         (GEOM_PATH_COL_KEY, fpath))

    logger.info("Writing species to {:s}".format(spc_csv_out))
    write_table_to_csv(spc_df, spc_csv_out)


def reactions_find_arrhenius(rxn_csv, rxn_csv_out, logger):
    """ get arrhenius parameters from job directories
    """
    from .ptorsscan import arrhenius as arrhenius_from_plog

    logger.info("Reading in {:s}".format(rxn_csv))
    rxn_df = pandas.read_csv(rxn_csv)

    col_keys = table_column_keys(rxn_df)
    assert 'reaction_id' in col_keys and 'path' in col_keys

    prefix = os.path.dirname(rxn_csv)

    def _get(rxn_row):
        arrh = None
        rid = rxn_row['reaction_id']
        if not is_empty_table_value(rxn_row['path']):
            path = os.path.join(prefix, rxn_row['path'])
            logger.info("reaction {:s}".format(rid))
            plog_path = os.path.join(path, 'rate.plog')
            if os.path.isfile(plog_path):
                logger.info(plog_path)
                plog_str = read_file(plog_path)
                arrh = arrhenius_from_plog(plog_str)
                if arrh:
                    logger.info("A={:f}, b={:f}, Ea={:f}".format(*arrh))
            else:
                logger.info("No rate.plog file found")
        return arrh if arrh else (None, None, None)

    arrh_col_keys = ('arrh_a', 'arrh_b', 'arrh_e')
    rxn_df = update_table_column_keys(rxn_df, col_keys=arrh_col_keys)
    arrh_as, arrh_bs, arrh_es = zip(
        *map(_get, iterate_table_row_dicts(rxn_df)))
    rxn_df['arrh_a'] = arrh_as
    rxn_df['arrh_b'] = arrh_bs
    rxn_df['arrh_e'] = arrh_es

    logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
    write_table_to_csv(rxn_df, rxn_csv_out)


def reactions_plot_arrhenius(rxn_csv, rxn_csv_ref, rxn_csv_out, plot_dir,
                             extension, tmp_rng, lbl_col_keys, id2path,
                             logger):
    """ make Arrhenius plots
    """
    from .plot import write_diagram
    from .iohelp import arrhenius_diagram

    logger.info("Reading in {:s}".format(rxn_csv))
    rxn_df = pandas.read_csv(rxn_csv)
    assert table_has_column_keys(rxn_df, ARRH_COL_KEYS)
    assert table_has_column_keys(rxn_df, lbl_col_keys)
    rxn_df = update_table_column_keys(rxn_df, ('plot_path',))

    assert len(tmp_rng) == 2 and tmp_rng[0] < tmp_rng[1]
    tmp_lo, tmp_hi = tmp_rng

    if rxn_csv_ref:
        logger.info("Reading in {:s}".format(rxn_csv_ref))
        rxn_df_ref = pandas.read_csv(rxn_csv_ref)
        assert table_has_column_keys(rxn_df_ref, ARRH_COL_KEYS)
    else:
        rxn_df_ref = None

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    for row in iterate_table_row_dicts(rxn_df):
        rid = row['reaction_id']
        logger.info("reaction {:s}".format(rid))

        cfts = tuple(map(row.__getitem__, ARRH_COL_KEYS))
        lbls = tuple(map(row.__getitem__, lbl_col_keys))
        ref_cfts = None
        if rxn_df_ref is not None:
            ref_row = table_lookup_row(rxn_df_ref, ('reaction_id', rid))
            if ref_row:
                ref_cfts = tuple(map(ref_row.__getitem__, ARRH_COL_KEYS))
        arrh_dgm = arrhenius_diagram(cfts, ref_cfts, tmp_lo, tmp_hi, lbls)

        if arrh_dgm:
            fname = '{:s}.{:s}'.format(id2path(rid), extension)
            fpath = os.path.join(plot_dir, fname)

            logger.info("  writing plot to {:s}".format(fpath))
            write_diagram(arrh_dgm, fpath, close=True)

            rxn_df = table_lookup_update(rxn_df,
                                         ('reaction_id', rid),
                                         ('plot_path', fpath))
        else:
            logger.info("  missing Arrhenius coefficients; skipping...")

    logger.info("Writing updated reaction table to {:s}".format(rxn_csv_out))
    write_table_to_csv(rxn_df, rxn_csv_out)


def read_thermo_data(spc_csv):
    """ a dictionary of thermo values (H298), indexed by species ID
    """
    thv_dct = None

    spc_df = pandas.read_csv(spc_csv)
    if 'therm_val' in spc_df:
        thv_dct = dict(zip(spc_df['species_id'], spc_df['therm_val']))

    return thv_dct


def chemkin_to_csv(mech_txt, thm_txt, rxn_csv_out, spc_csv_out, logger):
    """ parse CHEMKIN files
    """
    from .pchemkin import species as chemkin_species
    from .pchemkin import reactions as chemkin_reactions
    from .pchemkin import (thermodynamics_dictionaries as
                           chemkin_thermodynamics_dictionaries)
    from .pchemkin import kinetics as chemkin_kinetics

    logger.info("Reading in {:s}".format(mech_txt))
    mech_str = read_file(mech_txt)

    if thm_txt:
        logger.info("Reading in {:s}".format(thm_txt))
        thm_str = read_file(thm_txt)
    else:
        logger.info("No thermo file. Looking for thermo data in {:s}."
                    .format(mech_txt))
        thm_str = mech_str

    logger.info("Finding species")
    spcs = chemkin_species(mech_str)
    spc_ck_idxs = tuple(range(1, len(spcs)+1))

    logger.info("Finding reactions")
    rxns = chemkin_reactions(mech_str)
    rxn_ck_idxs = tuple(range(1, len(rxns)+1))

    logger.info("Finding thermodynamics data")
    thm_dcts = chemkin_thermodynamics_dictionaries(thm_str)

    logger.info("Finding kinetics data")
    kin_lst, reacs = chemkin_kinetics(mech_str)

    spc_df = table_from_columns((spc_ck_idxs, spcs),
                                ('chemkin_index', 'species'))
    rxn_df = table_from_columns((rxn_ck_idxs, rxns),
                                ('chemkin_index', 'reaction'))

    for rxn in rxns:
        if rxn not in reacs:
            logger.info(rxn)

    if kin_lst:
        assert len(kin_lst) == len(rxns)
        arrh_cols = tuple(zip(*kin_lst))
        rxn_df = append_table_columns(rxn_df, arrh_cols, ARRH_COL_KEYS)

    if thm_dcts:
        nasa_lo_dct, nasa_hi_dct, nasa_t_dct = thm_dcts
        nasa_lo_cols = tuple(zip(*map(nasa_lo_dct.__getitem__, spcs)))
        nasa_hi_cols = tuple(zip(*map(nasa_hi_dct.__getitem__, spcs)))
        nasa_t_cols = tuple(zip(*map(nasa_t_dct.__getitem__, spcs)))

        thm_col_keys = NASA_LO_COL_KEYS + NASA_HI_COL_KEYS + NASA_T_COL_KEYS
        thm_cols = nasa_lo_cols + nasa_hi_cols + nasa_t_cols
        spc_df = append_table_columns(spc_df, thm_cols, thm_col_keys)

    logger.info("Writing species to {:s}".format(spc_csv_out))
    write_table_to_csv(spc_df, spc_csv_out)

    logger.info("Writing reactions to {:s}".format(rxn_csv_out))
    write_table_to_csv(rxn_df, rxn_csv_out)


def chemkin_id_reactions(rxn_csv, spc_csv, rxn_csv_out, spc_csv_out, logger):
    """ determine reaction identifiers for CHEMKIN reactions
    """
    from .iohelp import translate_chemkin_reaction

    logger.info("Reading in {:s}".format(rxn_csv))
    rxn_df = pandas.read_csv(rxn_csv)

    logger.info("Reading in {:s}".format(spc_csv))
    spc_df = pandas.read_csv(spc_csv)

    assert table_has_column_keys(spc_df, ('species', 'species_id'))

    logger.info("Canonicalizing species IDs")
    sids = table_column(spc_df, 'species_id')
    can_sids = tuple(map(canonical_species_identifier, sids))
    spc_df = set_table_column(spc_df, 'species_id', can_sids)

    sid_dct = dict(zip(*table_columns(spc_df, ('species', 'species_id'))))

    rxns = table_column(rxn_df, 'reaction')
    rids = tuple(translate_chemkin_reaction(rxn, sid_dct) for rxn in rxns)
    can_rids = tuple(canonical_reaction_identifier(rid) if rid else None
                     for rid in rids)
    rxn_df = set_table_column(rxn_df, 'reaction_id', can_rids)

    spc_df = move_table_column_to_front(spc_df, 'species_id')
    rxn_df = move_table_column_to_front(rxn_df, 'reaction_id')

    logger.info("Writing species to {:s}".format(spc_csv_out))
    write_table_to_csv(spc_df, spc_csv_out)

    logger.info("Writing reactions to {:s}".format(rxn_csv_out))
    write_table_to_csv(rxn_df, rxn_csv_out)


def reactions_to_chemkin(cls, rxn_csv, spc_csv, mech_txt_out, logger):
    """ generate CHEMKIN files from CSVs
    """
    rct_sid_col_keys = dict(REACTANT_SID_COL_KEYS)[cls]
    prd_sid_col_keys = dict(PRODUCT_SID_COL_KEYS)[cls]

    logger.info("Reading in {:s}".format(spc_csv))
    spc_df = pandas.read_csv(spc_csv)
    spc_col_keys = table_column_keys(spc_df)
    assert 'species' in spc_col_keys and 'species_id' in spc_col_keys

    spc_dct = dict(zip(*table_columns(spc_df, ('species_id', 'species'))))

    def _chemkin_reaction_name(rct_sids, prd_sids):
        rct_str = '+'.join(map(spc_dct.__getitem__, rct_sids))
        prd_str = '+'.join(map(spc_dct.__getitem__, prd_sids))
        rxn = '='.join([rct_str, prd_str])
        return rxn

    logger.info("Reading in {:s}".format(rxn_csv))
    rxn_df = pandas.read_csv(rxn_csv)
    rxn_col_keys = table_column_keys(rxn_df)

    assert all(col_key in rxn_col_keys for col_key in rct_sid_col_keys)
    assert all(col_key in rxn_col_keys for col_key in prd_sid_col_keys)
    rct_sids_lst = tuple(zip(*table_columns(rxn_df, rct_sid_col_keys)))
    prd_sids_lst = tuple(zip(*table_columns(rxn_df, prd_sid_col_keys)))

    rxns = tuple(starmap(_chemkin_reaction_name,
                         zip(rct_sids_lst, prd_sids_lst)))

    assert all(col_key in rxn_col_keys for col_key in ARRH_COL_KEYS)
    arrh_cfts_lst = zip(*table_columns(rxn_df, ARRH_COL_KEYS))

    rxn_fmt = '{:{width}s} {:10.3e} {:8.3f} {:12.3f}'
    rxn_wd = max(map(len, rxns)) + 5
    format_ = partial(rxn_fmt.format, width=rxn_wd)
    rxn_block_str = '\n'.join(
        format_(rxn, *arrh_cfts) for rxn, arrh_cfts in zip(rxns, arrh_cfts_lst)
        if not any(map(is_empty_table_value, arrh_cfts)))

    mech_str = '\n'.join(['REACTIONS', rxn_block_str, 'END'])

    logger.info("Writing reactions to {:s}".format(mech_txt_out))
    write_file(mech_txt_out, mech_str)


def reactions_init(cls, rxn_csv, spc_csv, rxn_csv_out, cdt_csv_out, logger):
    """ initialize reactions
    """
    _init = reactions_initializer(
        cls=cls,
        is_candidate=dict(REACTION_CANDIDATE_FINDERS)[cls],
        reaction=dict(REACTION_FINDERS)[cls],
        sid_cols=dict(REACTION_SID_COL_KEYS)[cls],
        idx_cols=dict(REACTION_IDX_COL_KEYS)[cls]
    )
    return _init(spc_csv, rxn_csv, rxn_csv_out, cdt_csv_out, logger)


def reactions_run_old(cls, rxn_csv, spc_csv, rxn_rng_strs, tpl_txt, nodes,
                      run_dir, id2path, job_argv, logger):
    """ reactions runner
    """
    _run = reactions_runner(
        cls=cls,
        reaction_xyz_strings=dict(REACTION_XYZ_STRING_MAKERS)[cls],
        reaction_input_string=dict(REACTION_INPUT_STRING_MAKERS)[cls],
        sid_cols=dict(REACTION_SID_COL_KEYS)[cls],
        idx_cols=dict(REACTION_IDX_COL_KEYS)[cls]
    )
    _run(spc_csv, rxn_csv, tpl_txt, rxn_rng_strs, nodes, run_dir, id2path,
         job_argv, logger)


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
    rxn_df[key] = tuple(map(meets_condition_, rxn_df['reaction_id']))
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


def csv_intersect(table_csvs, col_key, table_csv_out, logger):
    """ intersect tables by column
    """
    table_dfs = []
    for table_csv in table_csvs:
        logger.info("Reading in {:s}".format(table_csv))
        table_df = pandas.read_csv(table_csv)
        table_dfs.append(table_df)

    table_df_out = intersect_tables(table_dfs, col_key)

    logger.info("Writing {:s}".format(table_csv_out))
    write_table_to_csv(table_df_out, table_csv_out)


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

        rxn_df_out = table_with_columns_like(rxn_df)
        cdt_df_out = table_with_columns_like(rxn_df)
        rxn_df_out = update_table_column_keys(rxn_df_out,
                                              col_keys=sid_cols+idx_cols)
        cdt_df_out = update_table_column_keys(cdt_df_out,
                                              col_keys=('exception',))

        rxn_df = update_table_column_keys(rxn_df, col_keys=('class',))

        for rxn_row in iterate_table_row_dicts(rxn_df):
            rid = rxn_row['reaction_id']
            if is_candidate(rid):
                logger.info('reaction: {:s}'.format(rid))

                err = None
                try:
                    rxn = reaction(rid, mgeo_dct, thv_dct)
                except Exception as err:
                    logger.info('  exception: {:s}!'.format(str(err)))
                    rxn = None

                if rxn:
                    sids, idxs = rxn

                    logger.info('  found {:s}!'.format(cls))

                    log_sids = ', '.join(
                        '{:s}: {:s}'.format(sid_col, sid)
                        for sid_col, sid in zip(sid_cols, sids))
                    log_idxs = ', '.join(
                        '{:s}: {:d}'.format(idx_col, idx)
                        for idx_col, idx in zip(idx_cols, idxs))

                    logger.info('  {:s}\n  {:s}'.format(log_sids, log_idxs))

                    rxn_df = table_lookup_update(rxn_df, ('reaction_id', rid),
                                                 ('class', cls))

                    rxn_row.update(zip(sid_cols, sids))
                    rxn_row.update(zip(idx_cols, idxs))
                    rxn_df_out = append_table_row_dicts(rxn_df_out, (rxn_row,))
                else:
                    rxn_row['exception'] = err
                    cdt_df_out = append_table_row_dicts(cdt_df_out, (rxn_row,))

        logger.info("Writing {:s} reactions to {:s}"
                    .format(cls, os.path.abspath(rxn_csv_out)))
        rxn_df_out = reindex_table(rxn_df_out)
        write_table_to_csv(rxn_df_out, rxn_csv_out)

        logger.info("Writing left-over candidates to {:s}"
                    .format(os.path.abspath(cdt_csv_out)))
        write_table_to_csv(cdt_df_out, cdt_csv_out)

        logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
        write_table_to_csv(rxn_df, rxn_csv)

    return _init


def reactions_setup_run(cls, rxn_csv, spc_csv, rxn_rng_strs, run_dir, id2path,
                        cmd_argv, logger):
    """ write xyz files for the runner
    """
    assert cls in ('abstraction', 'addition', 'migration')

    reaction_xyz_strings = dict(REACTION_XYZ_STRING_MAKERS)[cls]
    sid_cols = dict(REACTION_SID_COL_KEYS)[cls]
    idx_cols = dict(REACTION_IDX_COL_KEYS)[cls]

    logger.info("Reading in {:s}".format(rxn_csv))
    rxn_df = pandas.read_csv(rxn_csv)

    col_keys = table_column_keys(rxn_df)
    assert RID_COL_KEY in col_keys and RXN_IDX_COL_KEY in col_keys

    logger.info("Reading in species geometries from {:s}".format(spc_csv))
    mgeo_dct = read_geometries(spc_csv)

    if rxn_rng_strs:
        rxn_idxs = _interpret_range_strings(rxn_rng_strs)
        logger.info("Interpreted reaction index range argument: {:s}"
                    .format(repr(rxn_idxs)))
    else:
        logger.info("No reaction range argument. Running all reactions.")
        rxn_idxs = table_column(rxn_df, RXN_IDX_COL_KEY)

    if not os.path.exists(run_dir):
        logger.info("Creating run directory {:s}".format(run_dir))
        os.mkdir(run_dir)

    if not os.path.exists(run_dir):
        logger.info("Creating run directory {:s}".format(run_dir))
        os.mkdir(run_dir)

    def _create_job_dir(idx, rid, sids, idxs):
        logger.info("reaction {:d}: {:s}".format(idx, rid))
        logger.info('  indices: {:s}'.format(str(idxs)))

        dname = id2path(rid)
        dpath = os.path.join(run_dir, dname)
        logger.info("  Creating job directory {:s}".format(dpath))
        if not os.path.exists(dpath):
            os.mkdir(dpath)

        ret = (EMPTY, RXN_NOT_CREATED_VAL)

        dxyz_dct = reaction_xyz_strings(sids, idxs, mgeo_dct)
        if dxyz_dct:
            dxyz_sids = dxyz_dct.keys()
            dxyzs = dxyz_dct.values()

            fnames = tuple(map('{:s}.xyz'.format, map(id2path, dxyz_sids)))
            fpaths = tuple(os.path.join(dpath, fname) for fname in fnames)
            for fpath, dxyz in zip(fpaths, dxyzs):
                logger.info("  Writing {:s}".format(fpath))
                write_file(fpath, dxyz)

            ret = (dpath, RXN_CREATED_VAL)
        else:
            logger.info("  Failed to write .xyz files.")

        if cmd_argv:
            cmd_str = ' '.join(cmd_argv)
            logger.info("  running command `{:s}` in {:s}"
                        .format(cmd_str, dpath))
            try:
                subprocess.check_call(cmd_argv, cwd=dpath)
            except Exception as err:
                logger.info("  {:s}".format(err))

            logger.info('')

        return ret

    logger.info("Writing job .xyz files")

    sub_rxn_df = sql_where_in(rxn_df, RXN_IDX_COL_KEY, rxn_idxs)

    sub_idxs = tuple(sql_select_one(sub_rxn_df, RXN_IDX_COL_KEY))
    sub_rids = tuple(sql_select_one(sub_rxn_df, RID_COL_KEY))
    sub_sids_lst = tuple(zip(*(
        sql_select_one(sub_rxn_df, sid_col) for sid_col in sid_cols)))
    sub_idxs_lst = tuple(zip(*(
        sql_select_one(sub_rxn_df, idx_col) for idx_col in idx_cols)))

    paths, stats = zip(*starmap(
        _create_job_dir, zip(sub_idxs, sub_rids, sub_sids_lst, sub_idxs_lst)))

    rxn_df = update_table_column_keys(rxn_df,
                                      (RXN_PATH_COL_KEY, RXN_STAT_COL_KEY))

    rxn_df = update_column_by_index(rxn_df, row_indices(sub_rxn_df),
                                    RXN_PATH_COL_KEY, paths)
    rxn_df = update_column_by_index(rxn_df, row_indices(sub_rxn_df),
                                    RXN_STAT_COL_KEY, stats)

    logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
    write_table_to_csv(rxn_df, rxn_csv)


def reactions_run(cls, rxn_csv, rxn_rng_strs, tpl_txt, nodes, job_argv,
                  logger):
    """ reactions parallel runner
    """
    assert cls in ('abstraction', 'addition', 'migration')

    if not hasattr(job_argv, '__iter__'):
        raise ValueError("Missing run command.")

    sid_cols = dict(REACTION_SID_COL_KEYS)[cls]

    logger.info("Reading in {:s}".format(rxn_csv))
    rxn_df = pandas.read_csv(rxn_csv)

    col_keys = table_column_keys(rxn_df)
    assert (RID_COL_KEY in col_keys and RXN_IDX_COL_KEY in col_keys and
            RXN_PATH_COL_KEY in col_keys and RXN_STAT_COL_KEY in col_keys)

    if rxn_rng_strs:
        rxn_idxs = _interpret_range_strings(rxn_rng_strs)
        logger.info("Interpreted reaction index range argument: {:s}"
                    .format(repr(rxn_idxs)))
    else:
        logger.info("No reaction range argument. Running all reactions.")
        rxn_idxs = table_column(rxn_df, RXN_IDX_COL_KEY)

    logger.info("Reading template file from {:s}".format(tpl_txt))
    tpl_str = read_file(tpl_txt)

    def _submit_job(idx, rid, path, sids, worker_id):
        node = worker_id
        logger.info("reaction {:d}: {:s} assigned to node {:s} worker"
                    .format(idx, rid, node))

        inp_str = tpl_str.format(nodes=node, **dict(zip(sid_cols, sids)))
        inp_fpath = os.path.join(path, 'input.dat')

        logger.info("  node {:s} worker writing {:s}".format(node, inp_fpath))
        write_file(inp_fpath, inp_str)

        cmd_str = ' '.join(job_argv)
        logger.info("  node {:s} worker running {:s} in {:s}"
                    .format(node, cmd_str, path))

        ret = RXN_RAN_VAL
        try:
            subprocess.check_call(job_argv, cwd=path)
        except Exception as err:
            logger.info("  failed on node {:s} with error '{:s}'"
                        .format(node, err))
            ret = RXN_FAILED_VAL

        return ret

    run_rxn_df = sql_where_eq(
        sql_where_in(rxn_df, RXN_IDX_COL_KEY, rxn_idxs),
        RXN_STAT_COL_KEY,
        RXN_CREATED_VAL)

    run_idxs = tuple(sql_select_one(run_rxn_df, RXN_IDX_COL_KEY))
    run_rids = tuple(sql_select_one(run_rxn_df, RID_COL_KEY))
    run_paths = tuple(sql_select_one(run_rxn_df, RXN_PATH_COL_KEY))
    run_sids_lst = tuple(zip(*(
        sql_select_one(run_rxn_df, sid_col) for sid_col in sid_cols)))

    logger.info("Stepping through run directories to submit jobs")
    stats = tag_team_starmap(_submit_job,
                             zip(run_idxs, run_rids, run_paths, run_sids_lst),
                             worker_ids=nodes)
    rxn_df = update_column_by_index(rxn_df, row_indices(run_rxn_df),
                                    RXN_STAT_COL_KEY, stats)

    logger.info("Writing updated reaction table to {:s}".format(rxn_csv))
    write_table_to_csv(rxn_df, rxn_csv)


def reactions_runner(cls, reaction_xyz_strings, reaction_input_string,
                     sid_cols, idx_cols):
    """ run reactions
    """
    assert cls in ('abstraction', 'addition', 'migration')

    def _run(spc_csv, rxn_csv, tpl_txt, rxn_rng_strs, nodes, run_dir, id2path,
             job_argv, logger):
        logger.info("Reading in {:s}".format(rxn_csv))
        rxn_df = pandas.read_csv(rxn_csv)

        col_keys = table_column_keys(rxn_df)
        assert 'reaction_id' in col_keys and RXN_IDX_COL_KEY in col_keys

        logger.info("Reading in species geometries from {:s}".format(spc_csv))
        mgeo_dct = read_geometries(spc_csv)

        logger.info("Reading template file from {:s}".format(tpl_txt))
        tpl_str = read_file(tpl_txt)

        node_str = ', '.join(nodes)
        logger.info("Nodes: {:s}".format(node_str))
        tpl_keyval_dct = {'nodes': node_str}

        if rxn_rng_strs:
            rxn_idxs = _interpret_range_strings(rxn_rng_strs)
            logger.info("Interpreted reaction index range argument: {:s}"
                        .format(repr(rxn_idxs)))
        else:
            logger.info("No reaction range argument. Running all reactions.")
            rxn_idxs = table_column(rxn_df, RXN_IDX_COL_KEY)

        if not os.path.exists(run_dir):
            logger.info("Creating run directory {:s}".format(run_dir))
            os.mkdir(run_dir)

        logger.info("Writing job files")
        rxn_lkp = table_lookup_dictionary(rxn_df, RXN_IDX_COL_KEY)
        for idx in rxn_idxs:
            row = rxn_lkp[idx]
            rid = row['reaction_id']
            logger.info("reaction {:d}: {:s}".format(idx, rid))

            sids = tuple(map(row.__getitem__, sid_cols))
            idxs = tuple(map(row.__getitem__, idx_cols))
            logger.info('  indices: {:s}'.format(str(idxs)))

            dxyz_dct = reaction_xyz_strings(sids, idxs, mgeo_dct)
            if dxyz_dct:
                dxyz_sids = dxyz_dct.keys()
                dxyzs = dxyz_dct.values()

                dname = id2path(rid)
                dpath = os.path.join(run_dir, dname)
                logger.info("  Creating job directory {:s}".format(dpath))
                if not os.path.exists(dpath):
                    os.mkdir(dpath)

                fnames = tuple(map('{:s}.xyz'.format, map(id2path, dxyz_sids)))
                fpaths = tuple(os.path.join(dpath, fname) for fname in fnames)
                for fpath, dxyz in zip(fpaths, dxyzs):
                    logger.info("  Writing {:s}".format(fpath))
                    write_file(fpath, dxyz)

                inp_str = reaction_input_string(sids, tpl_str, tpl_keyval_dct)
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
        rxn_lkp = table_lookup_dictionary(rxn_df, RXN_IDX_COL_KEY)
        for idx in rxn_idxs:
            row = rxn_lkp[idx]
            if row['created']:
                rid = row['reaction_id']
                logger.info("reaction {:d}: {:s}".format(idx, rid))

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
    fle.write(str(contents))
    fle.close()


def read_file(file_txt):
    """ read file contents as a string
    """
    with open(file_txt, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()

    return file_str


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
        (strip_spaces(s) for s in split(':', kv))
        for kv in split('|', tmp_keyval_str))
    return tmp_keyval_dct


def _interpret_range_strings(rng_strs):

    def _interpret_range_string(rng_str):
        split_rng = split('-', rng_str)
        if len(split_rng) == 1:
            rng = [int(split_rng[-1])]
        elif len(split_rng) == 2:
            start, stop = map(int, split_rng)
            rng = list(range(start, stop+1))
        else:
            raise ValueError("Failed to interet index ranges")
        return rng

    return tuple(chain(*map(_interpret_range_string, rng_strs)))
