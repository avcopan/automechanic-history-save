""" drivers and functions for I/O
"""
import os
import json
import pandas


def write_file(fpath, contents, mode='w'):
    """ write contents to a file
    """
    fle = open(fpath, mode)
    fle.write(contents)
    fle.close()


def init_geometries(spc_csv, spc_csv_out, prefix, geom_dir, sid2fname,
                    logger):
    """ initialize species geometries
    """
    from .strid import xyz_string

    owd = os.getcwd()
    os.chdir(prefix)

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
    spc_df.to_csv(spc_csv_out)

    os.chdir(owd)


def init_from_rmg(rmg_mech_json, spc_csv_out, rxn_csv_out, prefix, logger):
    """ initialize the mechanism from RMG's JSON file
    """
    from .prmg import mechanism_species_identifiers
    from .prmg import mechanism_reaction_identifiers
    from .prmg import mechanism_uncertainties
    from .prmg import mechanism_sensitivities
    from .prmg import mechanism_importance_values

    logger.info("Parsing RMG mechanism JSON file")

    if not os.path.exists(prefix):
        os.mkdir(prefix)

    spc_csv_out_path = os.path.join(prefix, spc_csv_out)
    rxn_csv_out_path = os.path.join(prefix, rxn_csv_out)

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

    logger.info("Writing species to {:s}".format(spc_csv_out_path))
    spc_df.to_csv(spc_csv_out_path)
    logger.info("Writing reactions to {:s}".format(rxn_csv_out_path))
    rxn_df.to_csv(rxn_csv_out_path)


def init_from_chemkin(chemkin_mech_txt, spc_csv, spc_csv_out, rxn_csv_out,
                      prefix, logger):
    """ initialize the mechanism from a CHEMKIN file
    """
    from .pchemkin import mechanism_reaction_identifiers

    logger.info("Parsing CHEMKIN mechanism text file")

    if not os.path.exists(prefix):
        os.mkdir(prefix)

    spc_csv_out_path = os.path.join(prefix, spc_csv_out)
    rxn_csv_out_path = os.path.join(prefix, rxn_csv_out)

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

    logger.info("Writing species to {:s}".format(spc_csv_out_path))
    spc_df.to_csv(spc_csv_out_path)
    logger.info("Writing reactions to {:s}".format(rxn_csv_out_path))
    rxn_df.to_csv(rxn_csv_out_path)


def init_abstractions(spc_csv, rxn_csv):
    """ initialize hydrogen abstractions
    """
    pass


def init_additions(spc_csv, rxn_csv):
    """ initialize radical additions
    """
    pass
