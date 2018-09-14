""" functions for processing RMG data
"""
from .strid import reaction_identifier as reaction_identifier_from_sids


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


def species_thermo_value(spc_dct):
    """ species enthalpy at 298
    """
    return spc_dct['H298']


def reaction_name(rxn_dct):
    """ reaction name
    """
    return rxn_dct['name']


def reaction_identifier(rxn_dct):
    """ reaction ID
    """
    rct_sids, prd_sids = reaction_species_identifiers(rxn_dct)
    return reaction_identifier_from_sids(rct_sids, prd_sids)


def reaction_sensitivity(rxn_dct):
    """ reaction sensitivity
    """
    return rxn_dct['Sensitivity']


def reaction_value(rxn_dct):
    """ reaction value
    """
    return rxn_dct['Value']


def reaction_species_dictionaries(rxn_dct):
    """ return the species dictionaries for a reaction
    """
    return rxn_dct['Reactants'], rxn_dct['Products']


def reaction_species_identifiers(rxn_dct):
    """ species IDs for reactands and products
    """
    rct_spc_dcts, prd_spc_dcts = reaction_species_dictionaries(rxn_dct)
    rct_sids = tuple(map(species_identifier, rct_spc_dcts))
    prd_sids = tuple(map(species_identifier, prd_spc_dcts))
    return rct_sids, prd_sids
