""" coordinate library
"""
from .rot import unit_direction
from .rot import aligning_rotation_matrix
from .rot import local_coordinate_interpreter

NONSTEREO_XYZS = ((-1, 0, 0),         # atom 0 neighbor 0
                  (1, 0, 0),          # atom 0 neighbor 1
                  (0, 1, 0),          # atom 0 neighbor 2
                  (0, -1, 0))         # atom 0 neighbor 3

ATOM_STEREO_NEG_XYZS = ((-1, 0, 0),   # atom 0 neighbor 0
                        (0, 1, 0),    # atom 0 neighbor 1
                        (0, 0, -1),   # atom 0 neighbor 2
                        (0, -1, 0))   # atom 0 neighbor 3

ATOM_STEREO_POS_XYZS = ((-1, 0, 0),   # atom 0 neighbor 0
                        (0, 1, 0),    # atom 0 neighbor 1
                        (0, 0, +1),   # atom 0 neighbor 2
                        (0, -1, 0))   # atom 0 neighbor 3

BOND_STEREO_NEG_XYZS = ((-1, 0, 0),   # atom 0 neighbor 0
                        (1, 0, 0),    # atom 0 neighbor 1
                        (0, 1, 0),    # atom 1
                        (-1, 1, 0),   # atom 1 neighbor 0
                        (1, 1, 0))    # atom 1 neighbor 1

BOND_STEREO_POS_XYZS = ((-1, 0, 0),   # atom 0 neighbor 0
                        (1, 0, 0),    # atom 0 neighbor 1
                        (0, 1, 0),    # atom 1
                        (1, 1, 0),    # atom 1 neighbor 0
                        (-1, 1, 0))   # atom 1 neighbor 1


def nonstereo_coordinates(atm_xyz, anchor_xyz, anchor_key, sorted_keys):
    """ assign non-stereo coordinates from stencil
    """
    return coordinates_from_stencil(
        atm_xyz, anchor_xyz, anchor_key, sorted_keys,
        local_xyz_stencil=NONSTEREO_XYZS
    )


def atom_stereo_coordinates(atm_xyz, anchor_xyz, anchor_key, sorted_keys,
                            parity):
    """ assign atom-stereo coordinates from stencil
    """
    assert parity in (True, False)
    local_xyz_stencil = (ATOM_STEREO_POS_XYZS if parity else
                         ATOM_STEREO_NEG_XYZS)
    return coordinates_from_stencil(
        atm_xyz, anchor_xyz, anchor_key, sorted_keys,
        local_xyz_stencil=local_xyz_stencil
    )


def bond_stereo_coordinates(atm_xyz, anchor_xyz, anchor_key, sorted_keys,
                            parity):
    """ assign bond-stereo coordinates from stencil
    """
    assert parity in (True, False)
    local_xyz_stencil = (BOND_STEREO_POS_XYZS if parity else
                         BOND_STEREO_NEG_XYZS)
    return coordinates_from_stencil(
        atm_xyz, anchor_xyz, anchor_key, sorted_keys,
        local_xyz_stencil=local_xyz_stencil
    )


def coordinates_from_stencil(atm_xyz, anchor_xyz, anchor_key, sorted_keys,
                             local_xyz_stencil):
    """ coordinates for non-stereo atom neighbors
    """
    sorted_keys = list(sorted_keys)
    assert anchor_key in sorted_keys
    anchor_pos = sorted_keys.index(anchor_key)

    local_xyzs = list(local_xyz_stencil)

    anchor_key = sorted_keys.pop(anchor_pos)
    local_anchor_xyz = local_xyzs.pop(anchor_pos)

    local2global_ = local_coordinate_interpreter(
        trans=atm_xyz,
        rot=aligning_rotation_matrix(local_anchor_xyz,
                                     unit_direction(atm_xyz, anchor_xyz))
    )

    ngb_xyzs = list(map(local2global_, local_xyzs))

    assert len(sorted_keys) <= len(ngb_xyzs)
    return dict(zip(sorted_keys, ngb_xyzs))
