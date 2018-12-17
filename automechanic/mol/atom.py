""" atom functions
"""


def valence(atm):
    """ bonding valence
    """
    return {'H': 1,
            'C': 4,
            'N': 3,
            'O': 2, 'S': 2,
            'F': 1, 'CL': 1,
            'HE': 0, 'AR': 0}[atm.upper()]


def lone_pair_count(atm):
    """ lone pair count
    """
    return {'H': 0,
            'C': 0,
            'N': 1,
            'O': 2, 'S': 2,
            'F': 3, 'CL': 3,
            'HE': 1, 'AR': 4}[atm.upper()]
