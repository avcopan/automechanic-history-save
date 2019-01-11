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
