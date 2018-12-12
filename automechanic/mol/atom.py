""" atom functions
"""


def valence(atm):
    """ valence of an atom
    """
    vlnc_dct = {'H': 1,
                'C': 4,
                'N': 3,
                'O': 2, 'S': 2,
                'F': 1, 'Cl': 1,
                'HE': 0, 'AR': 0}
    return vlnc_dct[atm.upper()]
