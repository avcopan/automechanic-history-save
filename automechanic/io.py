""" messy I/O based drivers and functions
"""
import os
import sys
import subprocess
from future.moves.itertools import filterfalse
import pandas
from . import smiles
from . import dotxyz
from . import pchemkin
from . import geomlib
from . import smileslib
from .timeout import TimeoutError


def write_file(fpath, contents, mode='w'):
    """ write contents to a file
    """
    fle = open(fpath, mode)
    fle.write(contents)
    fle.close()


def write_geometries(spc_df, geom_dir, smi2fname):
    """ write geometries to .xyz files
    """
    assert isinstance(spc_df, pandas.DataFrame) and 'smiles' in spc_df

    if 'path' not in spc_df:
        spc_df['path'] = None

    if not os.path.exists(geom_dir):
        os.makedirs(geom_dir)

    for idx, row in spc_df.iterrows():
        smi = row['smiles']
        fname = smi2fname(smi)
        fpath = os.path.join(geom_dir, fname)
        if not os.path.exists(fpath):
            write_file(fpath, contents=smiles.xyz_string(smi))
        if spc_df.at[idx, 'path'] is None:
            spc_df.at[idx, 'path'] = fpath

    return spc_df


def read_geometries(spc_df):
    """ read geometries from .xyz files
    """
    assert (isinstance(spc_df, pandas.DataFrame) and 'smiles' in spc_df
            and 'path' in spc_df)

    mgeo_dct = {}

    for idx, row in spc_df.iterrows():
        smi = row['smiles']
        fpath = row['path']
        assert os.path.exists(fpath)
        dxyz = open(fpath).read()
        mgeo = dotxyz.geometry(dxyz)
        mgeo_dct[smi] = mgeo

    return mgeo_dct


def parse_mechanism(mech_str, spc_df, excl=('OHV', 'CHV')):
    """ parse the mechanism in a CHEMKIN file
    """
    smi_dct = dict(zip(spc_df['species'], spc_df['smiles']))

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

    def _to_smirks(rxn):
        rspcs, pspcs = rxn
        rsmis = [smi_dct[spc] for spc in rspcs]
        psmis = [smi_dct[spc] for spc in pspcs]
        smrk = smiles.make_smirks(rsmis, psmis)
        return smrk

    pi_smrks = list(map(_to_smirks, incl_pi_rxns))
    lp_smrks = list(map(_to_smirks, incl_lp_rxns))
    fo_smrks = list(map(_to_smirks, incl_fo_rxns))

    pi_df = pandas.DataFrame({'smirks': pi_smrks})
    lp_df = pandas.DataFrame({'smirks': lp_smrks})
    fo_df = pandas.DataFrame({'smirks': fo_smrks})
    rxn_df = pandas.concat([pi_df, lp_df, fo_df],
                           keys=['press_indep', 'low_p', 'falloff'])
    return rxn_df


def find_potential_hydrogen_abstractions(rxn_df):
    """ find potential hydrogen abstractions
    """
    print('finding potential hydrogen abstractions by formula...')
    assert isinstance(rxn_df, pandas.DataFrame) and 'smirks' in rxn_df

    if 'maybe_abstr' not in rxn_df:
        rxn_df['maybe_abstr'] = None

    rxn_df = rxn_df.rename(columns={"smirks": "orig_smirks"})
    rxn_df['smirks'] = None

    for idx, row in rxn_df.iterrows():
        orig_smrk = row['orig_smirks']
        sys.stdout.flush()
        smrk = smileslib.match_hydrogen_abstraction_formula(orig_smrk)
        if smrk:
            rxn_df.at[idx, 'maybe_abstr'] = True
        else:
            smrk = orig_smrk
        rxn_df.at[idx, 'smirks'] = smrk
    return rxn_df


def initialize_hydrogen_abstractions(rxn_df, spc_df, abstr_dir, smi2fname,
                                     smrk2dname):
    """ find hydrogen abstractions in the reaction database
    """
    print('reading geometries from .xyz files...')
    mgeo_dct = read_geometries(spc_df)

    if 'type' not in rxn_df:
        rxn_df['type'] = None
    if 'abstr_exception' not in rxn_df:
        rxn_df['abstr_exception'] = None

    if not os.path.exists(abstr_dir):
        os.makedirs(abstr_dir)

    abstr_df = pandas.DataFrame(columns=['smirks', 'path'])

    print('looping over potential hydrogen abstractions...')
    for idx, row in rxn_df.iterrows():
        if rxn_df.at[idx, 'maybe_abstr']:
            smrk = row['smirks']
            print('reaction {:d}: {:s}'.format(idx, smrk))
            rsmis, psmis = smiles.split_smirks(smrk)
            rmgeos = tuple(mgeo_dct[smi] for smi in rsmis)
            pmgeos = tuple(mgeo_dct[smi] for smi in psmis)
            try:
                abstr = geomlib.find_hydrogen_abstraction(rmgeos, pmgeos)
            except RuntimeError as err:
                rxn_df.at[idx, 'abstr_exception'] = 'RuntimeError'
                print('  abstraction finder failed: {:s}'.format(str(err)))
                continue
            except TimeoutError as err:
                rxn_df.at[idx, 'abstr_exception'] = 'TimeoutError'
                print('  abstraction finder failed: {:s}'.format(str(err)))
                continue
            if abstr:
                print('  found hydrogen abstraction...')
                rxn_df.at[idx, 'type'] = 'abstraction'
                forw, back = abstr
                (q1h_mgeo, q1h_idx), (q2_mgeo, q2_idx) = forw
                (q1_mgeo, _), (q2h_mgeo, _) = back

                q1h_smi, q2_smi = rsmis
                q1_smi, q2h_smi = psmis

                try:
                    q1h_dxyz = geomlib.abstractee_xyz_string(q1h_mgeo, q1h_idx)
                except ValueError as err:
                    rxn_df.at[idx, 'abstr_exception'] = 'ValueError'
                    print('  failed to generate abstractee xyz: {:s}'
                          .format(str(err)))
                    continue

                q2_dxyz = geomlib.abstractor_xyz_string(q2_mgeo, q2_idx)
                q1_dxyz = smiles.xyz_string(q1_smi)
                q2h_dxyz = smiles.xyz_string(q2h_smi)

                rxn_dname = smrk2dname(smrk)
                rxn_path = os.path.join(abstr_dir, rxn_dname)

                if not os.path.exists(rxn_path):
                    os.makedirs(rxn_path)

                q1h_fname = smi2fname(q1h_smi)
                q2_fname = smi2fname(q2_smi)
                q1_fname = smi2fname(q1_smi)
                q2h_fname = smi2fname(q2h_smi)

                print('  writing xyz files...')
                q1h_fpath = os.path.join(rxn_path, q1h_fname)
                q2_fpath = os.path.join(rxn_path, q2_fname)
                q1_fpath = os.path.join(rxn_path, q1_fname)
                q2h_fpath = os.path.join(rxn_path, q2h_fname)

                write_file(q1h_fpath, q1h_dxyz)
                write_file(q2_fpath, q2_dxyz)
                write_file(q1_fpath, q1_dxyz)
                write_file(q2h_fpath, q2h_dxyz)

                abstr_df = abstr_df.append({'smirks': smrk, 'path': rxn_path},
                                           ignore_index=True)

    return rxn_df, abstr_df
