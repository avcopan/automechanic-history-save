""" functions for parsing RMG files
"""
from itertools import starmap
from itertools import chain
from .strid import reaction_identifier


def species_name(spc_dct):
    """ species name
    """
    return spc_dct['name']


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


def mechanism_species_dictionaries(mech_rxn_dcts):
    """ dictionaries for species in the mechanism
    """
    names = []
    spc_dcts = []
    for rxn_dct in mech_rxn_dcts:
        for spc_dct in chain(rxn_dct['Reactants'], rxn_dct['Products']):
            name = species_name(spc_dct)
            if name not in names:
                names.append(name)
                spc_dcts.append(spc_dct)
    return spc_dcts
#
#
# def mechanism_species_identifiers(mech_rxn_dcts):
#     """ IDs for species in the mechanism
#     """
#     rct_sids_lst, prd_sids_lst = zip(*map(reaction_species_identifiers,
#                                           mech_rxn_dcts))
#     sids = reduce(set.union, map(set, chain(rct_sids_lst, prd_sids_lst)))
#     sids = sorted(sorted(sids), key=len)
#     return sids


def mechanism_reaction_identifiers(mech_rxn_dcts):
    """ IDs for reactions in the mechanism
    """
    rct_sids_lst, prd_sids_lst = zip(*map(reaction_species_identifiers,
                                          mech_rxn_dcts))
    rids = tuple(starmap(reaction_identifier, zip(rct_sids_lst, prd_sids_lst)))
    return rids


def mechanism_sensitivities(mech_rxn_dcts):
    """ reaction uncertainties
    """
    stvts = tuple(rxn_dct['Sensitivity'] for rxn_dct in mech_rxn_dcts)
    return stvts


def mechanism_importance_values(mech_rxn_dcts):
    """ reaction uncertainties
    """
    ipvls = tuple(rxn_dct['Value'] for rxn_dct in mech_rxn_dcts)
    return ipvls


def mechanism_reaction_names(mech_rxn_dcts):
    """ reactions names
    """
    names = tuple(rxn_dct['name'] for rxn_dct in mech_rxn_dcts)
    return names
