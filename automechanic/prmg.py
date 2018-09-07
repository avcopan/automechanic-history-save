""" functions for parsing RMG files
"""
from itertools import starmap
from itertools import chain
from functools import reduce
from .strid import reaction_identifier


def species_identifier(spc_dct):
    """ determine the species identifier from an RMG species dictionary
    """
    smi = spc_dct['SMILES'][0]
    mult = spc_dct['multiplicity']
    sid = '{:s}_m{:d}'.format(smi, mult)
    return sid


def reaction_species_identifiers(rxn_dct):
    """ parse an RMG reaction dictionary
    """
    rct_dcts = rxn_dct['Reactants']
    prd_dcts = rxn_dct['Products']

    rct_sids = list(map(species_identifier, rct_dcts))
    prd_sids = list(map(species_identifier, prd_dcts))

    return rct_sids, prd_sids


def mechanism_species_identifiers(mech_rxn_dcts):
    """ IDs for species in the mechanism
    """
    rct_sids_lst, prd_sids_lst = zip(*map(reaction_species_identifiers,
                                          mech_rxn_dcts))
    sids = reduce(set.union, map(set, chain(rct_sids_lst, prd_sids_lst)))
    sids = sorted(sorted(sids), key=len)
    return sids


def mechanism_reaction_identifiers(mech_rxn_dcts):
    """ IDs for reactions in the mechanism
    """
    rct_sids_lst, prd_sids_lst = zip(*map(reaction_species_identifiers,
                                          mech_rxn_dcts))
    rids = tuple(starmap(reaction_identifier, zip(rct_sids_lst, prd_sids_lst)))
    return rids


# def mechanism_uncertainties(mech_rxn_dcts):
#     """ reaction uncertainties
#     """
#     ucrts = tuple(rxn_dct['Uncertainty'] for rxn_dct in mech_rxn_dcts)
#     return ucrts
#
#
# def mechanism_sensitivities(mech_rxn_dcts):
#     """ reaction uncertainties
#     """
#     stvts = tuple(rxn_dct['Sensitivity'] for rxn_dct in mech_rxn_dcts)
#     return stvts


def mechanism_importance_values(mech_rxn_dcts):
    """ reaction uncertainties
    """
    ipvls = tuple(rxn_dct['Value'] for rxn_dct in mech_rxn_dcts)
    return ipvls
