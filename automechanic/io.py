""" messy I/O based drivers and functions
"""
import os
import sys
from future.moves.itertools import filterfalse
import pandas
from . import strid
from . import dotxyz
from . import pchemkin
from . import geom
from . import geomlib
from . import stridlib
from .timeout import TimeoutError


def write_file(fpath, contents, mode='w'):
    """ write contents to a file
    """
    fle = open(fpath, mode)
    fle.write(contents)
    fle.close()


def write_geometries(spc_df, path, sid2fname):
    """ write geometries to .xyz files
    """
    assert isinstance(spc_df, pandas.DataFrame) and 'species_id' in spc_df

    if 'path' not in spc_df:
        spc_df['path'] = None

    if not os.path.exists(path):
        os.makedirs(path)

    for idx, row in spc_df.iterrows():
        sid = row['species_id']
        fname = sid2fname(sid)
        fpath = os.path.join(path, fname)
        if not os.path.exists(fpath):
            dxyz = strid.xyz_string(sid)
            write_file(fpath, contents=dxyz)
        if spc_df.at[idx, 'path'] is None:
            spc_df.at[idx, 'path'] = fpath

    return spc_df


def read_geometries(spc_df):
    """ read geometries from .xyz files
    """
    assert (isinstance(spc_df, pandas.DataFrame) and 'species_id' in spc_df
            and 'path' in spc_df)

    mgeo_dct = {}

    for idx, row in spc_df.iterrows():
        sid = row['species_id']
        fpath = row['path']
        assert os.path.exists(fpath)
        dxyz = open(fpath).read()
        mgeo = dotxyz.geometry(dxyz)
        mgeo_dct[sid] = mgeo

    return mgeo_dct


def parse_mechanism(mech_str, spc_df, excl=('OHV', 'CHV')):
    """ parse the mechanism in a CHEMKIN file
    """
    sid_dct = dict(zip(spc_df['species'], spc_df['species_id']))

    def _contains_species(spc_lst):
        def __yes(rxn):
            rcts, prds = rxn
            rgts = rcts + prds
            return any(s in rgts for s in spc_lst) if spc_lst else False
        return __yes

    pi_rxns, lp_rxns, fo_rxns = pchemkin.find_reactions(mech_str)

    contains_excluded = _contains_species(excl)
    incl_pi_rxns = list(filterfalse(contains_excluded, pi_rxns))
    incl_lp_rxns = list(filterfalse(contains_excluded, lp_rxns))
    incl_fo_rxns = list(filterfalse(contains_excluded, fo_rxns))

    def _to_rid(rxn):
        rct_spcs, prd_spcs = rxn
        rct_sids = tuple(map(sid_dct.__getitem__, rct_spcs))
        prd_sids = tuple(map(sid_dct.__getitem__, prd_spcs))
        rid = strid.reaction_identifier(rct_sids, prd_sids)
        return rid

    pi_rids = list(map(_to_rid, incl_pi_rxns))
    lp_rids = list(map(_to_rid, incl_lp_rxns))
    fo_rids = list(map(_to_rid, incl_fo_rxns))

    pi_df = pandas.DataFrame({'reaction_id': pi_rids})
    lp_df = pandas.DataFrame({'reaction_id': lp_rids})
    fo_df = pandas.DataFrame({'reaction_id': fo_rids})
    rxn_df = pandas.concat([pi_df, lp_df, fo_df],
                           keys=['press_indep', 'low_p', 'falloff'],
                           names=['pressure dependence', ''])
    return rxn_df


def find_potential_hydrogen_abstractions(rxn_df):
    """ find potential hydrogen abstractions
    """
    print('finding potential hydrogen abstractions by formula...')
    assert isinstance(rxn_df, pandas.DataFrame) and 'reaction_id' in rxn_df

    if 'maybe_abstr' not in rxn_df:
        rxn_df['maybe_abstr'] = None

    rxn_df = rxn_df.rename(columns={"reaction_id": "orig_reaction_id"})
    rxn_df['reaction_id'] = None

    for idx, row in rxn_df.iterrows():
        orig_rid = row['orig_reaction_id']
        sys.stdout.flush()
        rid = stridlib.match_hydrogen_abstraction_formula(orig_rid)
        if rid:
            rxn_df.at[idx, 'maybe_abstr'] = True
        else:
            rid = orig_rid
            rxn_df.at[idx, 'maybe_abstr'] = False
        rxn_df.at[idx, 'reaction_id'] = rid
    return rxn_df


def initialize_hydrogen_abstractions(rxn_df, spc_df, path, sid2fname,
                                     rid2dname):
    """ find hydrogen abstractions in the reaction database
    """
    print('reading geometries from .xyz files...')
    mgeo_dct = read_geometries(spc_df)

    if 'type' not in rxn_df:
        rxn_df['type'] = None
    if 'abstr_exception' not in rxn_df:
        rxn_df['abstr_exception'] = None

    if not os.path.exists(path):
        os.makedirs(path)

    abstr_df = pandas.DataFrame(columns=['reaction_id', 'path'])

    print('looping over potential hydrogen abstractions...')
    for idx, row in rxn_df.iterrows():
        if rxn_df.at[idx, 'maybe_abstr']:
            rid = row['reaction_id']
            print('reaction {:d}: {:s}'.format(idx, rid))

            rct_sids, prd_sids = strid.split_reaction_identifier(rid)
            q1h_sid, q2_sid = rct_sids
            q1_sid, q2h_sid = prd_sids

            q1h_mgeo = mgeo_dct[q1h_sid]
            q2_mgeo = mgeo_dct[q2_sid]
            q1_mgeo = mgeo_dct[q1_sid]
            q2h_mgeo = mgeo_dct[q2h_sid]

            try:
                idxs = geomlib.find_hydrogen_abstraction(q1h_mgeo, q2_mgeo,
                                                         q1_mgeo, q2h_mgeo)
            except RuntimeError as err:
                rxn_df.at[idx, 'abstr_exception'] = 'RuntimeError'
                print('  abstraction finder failed: {:s}'.format(str(err)))
                continue
            except TimeoutError as err:
                rxn_df.at[idx, 'abstr_exception'] = 'TimeoutError'
                print('  abstraction finder failed: {:s}'.format(str(err)))
                continue

            if idxs:
                print('  found hydrogen abstraction...')
                rxn_df.at[idx, 'type'] = 'abstraction'

                q1h_idx, q2_idx, _, _ = idxs

                try:
                    q1h_dxyz = geomlib.abstractee_xyz_string(q1h_mgeo, q1h_idx)
                except ValueError as err:
                    rxn_df.at[idx, 'abstr_exception'] = 'ValueError'
                    print('  failed to generate abstractee xyz: {:s}'
                          .format(str(err)))
                    continue

                q2_dxyz = geomlib.abstractor_xyz_string(q2_mgeo, q2_idx)
                q1_dxyz = geom.xyz_string(q1_mgeo)
                q2h_dxyz = geom.xyz_string(q2h_mgeo)

                dname = rid2dname(rid)
                dpath = os.path.join(path, dname)

                if not os.path.exists(dpath):
                    os.makedirs(dpath)

                q1h_fname = sid2fname(q1h_sid)
                q2_fname = sid2fname(q2_sid)
                q1_fname = sid2fname(q1_sid)
                q2h_fname = sid2fname(q2h_sid)

                print('  writing xyz files...')
                q1h_fpath = os.path.join(dpath, q1h_fname)
                q2_fpath = os.path.join(dpath, q2_fname)
                q1_fpath = os.path.join(dpath, q1_fname)
                q2h_fpath = os.path.join(dpath, q2h_fname)

                write_file(q1h_fpath, q1h_dxyz)
                write_file(q2_fpath, q2_dxyz)
                write_file(q1_fpath, q1_dxyz)
                write_file(q2h_fpath, q2h_dxyz)

                abstr_df = abstr_df.append({'reaction_id': rid, 'path': dpath},
                                           ignore_index=True)

    return rxn_df, abstr_df
