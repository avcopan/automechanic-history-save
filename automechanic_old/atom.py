""" atom functions
"""


def valence(atm):
    """ valence of an atom
    """
    vlnc_dct = {'H': 1, 'C': 4, 'N': 3, 'O': 2, 'S': 2, 'Cl': 1}
    return vlnc_dct[atm]
