""" coordinate library
"""

NON_STEREO_ATOM_NEIGHBORS = ((-1, 0, 0),    # atom 0
                             (1, 0, 0),     # atom 1
                             (0, 1, 0),     # atom 2
                             (0, -1, 0))    # atom 3

NEG_STEREO_ATOM_NEIGHBORS = ((-1, 0, 0),     # atom 0
                             (0, 1, 0),      # atom 1
                             (0, 0, -1),     # atom 2
                             (0, -1, 0))     # atom 3

POS_STEREO_ATOM_NEIGHBORS = ((-1, 0, 0),     # atom 0
                             (0, 1, 0),      # atom 1
                             (0, 0, +1),     # atom 2
                             (0, -1, 0))     # atom 3


def assign_non_stereo_atom_neighbors(atm_key, ngb_key, ngb_keys, atm_xyz_dct,
                                     key_sorter):
    """ assign non-stereo atom neighbors
    """
    print(atm_key)
    print(ngb_key)
    print(ngb_keys)
    print(atm_xyz_dct)
    print(key_sorter)
