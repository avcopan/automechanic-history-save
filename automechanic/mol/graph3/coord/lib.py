""" coordinate library
"""
from itertools import chain as _chain
from .rot import unit_direction
from .rot import aligning_rotation_matrix
from .rot import local_coordinate_interpreter

NONSTEREO_XYZS = ((-1, 0, 0),         # atom 0 neighbor 0
                  (1, 0, 0),          # atom 0 neighbor 1
                  (0, 1, 0),          # atom 0 neighbor 2
                  (0, -1, 0))         # atom 0 neighbor 3


def nonstereo_coordinates(atm_xyz, anchor_xyz, anchor_key, sorted_keys):
    """ assign non-stereo coordinates from stencil
    """
    return coordinates_from_stencil(
        atm_xyz, anchor_xyz, anchor_key, sorted_keys,
        local_xyz_stencil=NONSTEREO_XYZS
    )


# def atom_stereo_coordinates(atm_xyz, anchor_xyz, anchor_key, sorted_keys,
#                             parity):
#     """ assign atom-stereo coordinates from stencil
#     """
#     assert parity in (True, False)
#     local_xyz_stencil = (ATOM_STEREO_POS_XYZS if parity else
#                          ATOM_STEREO_NEG_XYZS)
#     return coordinates_from_stencil(
#         atm_xyz, anchor_xyz, anchor_key, sorted_keys,
#         local_xyz_stencil=local_xyz_stencil
#     )


def atom_stereo_coordinates(anchor_key, atm_key, atm_ngb_keys_dct, xyz_dct,
                            key_sorter, parity):
    """ assign atom-stereo coordinates from stencil
    """
    local_xyz_stencil = ((0, 0, 0),                     # atm 1
                         (-1, 0, 0),                    # atm 1 ngb 0
                         (0, 1, 0),                     # atm 1 ngb 1
                         (0, 0, (-1)**(not parity)),    # atm 1 ngb 2
                         (0, -1, 0))                    # atm 1 ngb 3

    print(anchor_key)
    print(atm_key)
    print(atm_ngb_keys_dct)
    print(xyz_dct)
    print(key_sorter)
    print(parity)
    print(local_xyz_stencil)


def bond_stereo_coordinates(anchor_key, bnd_key, atm_ngb_keys_dct, xyz_dct,
                            key_sorter, parity):
    """ assign bond-stereo coordinates from stencil
    """
    local_xyz_stencil = ((0, 0, 0),                     # atm 1
                         (0, 1, 0),                     # atm 2
                         (-1, 0, 0),                    # atm 1 ngb 0 (anchor?)
                         (1, 0, 0),                     # atm 1 ngb 1 (anchor?)
                         ((-1) ** (not parity), 1, 0),  # atm 2 ngb 0
                         ((-1) ** parity, 1, 0))        # atm 2 ngb 1

    atm1_key, atm2_key = bnd_key
    atm1_ngb_keys = atm_ngb_keys_dct[atm1_key]
    atm2_ngb_keys = atm_ngb_keys_dct[atm2_key]

    if anchor_key in atm2_ngb_keys:
        atm1_key, atm2_key = atm2_key, atm1_key
        atm1_ngb_keys, atm2_ngb_keys = atm2_ngb_keys, atm1_ngb_keys

    sorted_keys = list(_chain(
        [atm1_key, atm2_key],
        key_sorter(filter(lambda x: x != atm2_key, atm1_ngb_keys)),
        key_sorter(filter(lambda x: x != atm1_key, atm2_ngb_keys))))

    assert atm1_key in xyz_dct
    assert anchor_key in xyz_dct
    xyz_dct = coordinates_from_stencil(
        atm_xyz=xyz_dct[atm1_key],
        anchor_xyz=xyz_dct[anchor_key],
        anchor_key=anchor_key,
        sorted_keys=sorted_keys,
        local_xyz_stencil=local_xyz_stencil
    )
    return xyz_dct


def coordinates_from_stencil(atm_xyz, anchor_xyz, anchor_key, sorted_keys,
                             local_xyz_stencil):
    """ coordinates for non-stereo atom neighbors
    """
    sorted_keys = list(sorted_keys)
    assert anchor_key in sorted_keys
    anchor_pos = sorted_keys.index(anchor_key)

    local_anchor_xyz = local_xyz_stencil[anchor_pos]

    local2global_ = local_coordinate_interpreter(
        trans=atm_xyz,
        rot=aligning_rotation_matrix(local_anchor_xyz,
                                     unit_direction(atm_xyz, anchor_xyz))
    )

    ngb_xyzs = list(map(local2global_, local_xyz_stencil))

    assert len(sorted_keys) <= len(ngb_xyzs)
    return dict(zip(sorted_keys, ngb_xyzs))
