""" coordinate library
"""
from itertools import chain as _chain
from .linalg import unit_direction
from .linalg import aligning_rotation_matrix
from .linalg import local_coordinate_interpreter


def nonstereo_coordinates(anchor_key, atm_key, atm_ngb_keys_dct, xyz_dct):
    """ assign non-stereo coordinates from a stencil
    """
    stencil_xyzs = ((0, 0, 0),    # atm 1
                    (-1, 0, 0),   # atm 1 ngb 0
                    (1, 0, 0),    # atm 1 ngb 1
                    (0, 1, 0),    # atm 1 ngb 2
                    (0, -1, 0))   # atm 1 ngb 3

    assert atm_key in atm_ngb_keys_dct
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    stencil_keys = list(_chain([atm_key], atm_ngb_keys))

    xyz_dct = dict.copy(xyz_dct)
    xyz_dct.update(_from_stencil(
        atm_key, anchor_key, xyz_dct, stencil_keys, stencil_xyzs))

    boundary_edges = tuple((atm_key, ngb_key) for ngb_key in atm_ngb_keys
                           if ngb_key != anchor_key)

    return xyz_dct, boundary_edges


def atom_stereo_coordinates(anchor_key, atm_key, atm_ngb_keys_dct, xyz_dct,
                            key_sorter, parity):
    """ assign atom-stereo coordinates from a stencil
    """
    stencil_xyzs = ((0, 0, 0),                     # atm 1
                    (-1, 0, 0),                    # atm 1 ngb 0
                    (0, 1, 0),                     # atm 1 ngb 1
                    (0, 0, (-1)**(not parity)),    # atm 1 ngb 2
                    (0, -1, 0))                    # atm 1 ngb 3

    assert atm_key in atm_ngb_keys_dct
    atm_ngb_keys = atm_ngb_keys_dct[atm_key]

    stencil_keys = list(_chain([atm_key], key_sorter(atm_ngb_keys)))

    assert len(stencil_keys) == len(stencil_xyzs)
    xyz_dct = dict.copy(xyz_dct)
    xyz_dct.update(_from_stencil(
        atm_key, anchor_key, xyz_dct, stencil_keys, stencil_xyzs))

    boundary_edges = tuple((atm_key, ngb_key) for ngb_key in atm_ngb_keys
                           if ngb_key != anchor_key)

    return xyz_dct, boundary_edges


def bond_stereo_coordinates(anchor_key, bnd_key, atm_ngb_keys_dct, xyz_dct,
                            key_sorter, parity):
    """ assign bond-stereo coordinates from a stencil
    """
    stencil_xyzs = ((0, 0, 0),                     # atm 1
                    (0, 1, 0),                     # atm 2
                    (-1, 0, 0),                    # atm 1 ngb 0 (anchor?)
                    (1, 0, 0),                     # atm 1 ngb 1 (anchor?)
                    ((-1) ** (not parity), 1, 0),  # atm 2 ngb 0
                    ((-1) ** parity, 1, 0))        # atm 2 ngb 1

    atm1_key, atm2_key = bnd_key
    assert atm1_key in atm_ngb_keys_dct
    assert atm2_key in atm_ngb_keys_dct
    atm1_ngb_keys = atm_ngb_keys_dct[atm1_key]
    atm2_ngb_keys = atm_ngb_keys_dct[atm2_key]

    if anchor_key in atm2_ngb_keys:
        atm1_key, atm2_key = atm2_key, atm1_key
        atm1_ngb_keys, atm2_ngb_keys = atm2_ngb_keys, atm1_ngb_keys

    stencil_keys = list(_chain(
        [atm1_key, atm2_key],
        key_sorter(filter(lambda x: x != atm2_key, atm1_ngb_keys)),
        key_sorter(filter(lambda x: x != atm1_key, atm2_ngb_keys))))

    assert len(stencil_keys) == len(stencil_xyzs)
    xyz_dct = dict.copy(xyz_dct)
    xyz_dct.update(_from_stencil(
        atm1_key, anchor_key, xyz_dct, stencil_keys, stencil_xyzs))

    boundary_edges = tuple(_chain(
        ((atm1_key, ngb_key) for ngb_key in atm1_ngb_keys
         if ngb_key not in (anchor_key, atm2_key)),
        ((atm2_key, ngb_key) for ngb_key in atm2_ngb_keys
         if ngb_key != atm1_key)))

    return xyz_dct, boundary_edges


def _from_stencil(atm_key, anchor_key, xyz_dct, stencil_keys, stencil_xyzs):
    assert atm_key in xyz_dct
    assert anchor_key in xyz_dct

    true_atm_xyz = xyz_dct[atm_key]
    true_anchor_xyz = xyz_dct[anchor_key]
    true_bond_xyz = unit_direction(true_atm_xyz, true_anchor_xyz)

    assert atm_key in stencil_keys
    assert anchor_key in stencil_keys
    stencil_atm_xyz = stencil_xyzs[stencil_keys.index(atm_key)]
    stencil_anchor_xyz = stencil_xyzs[stencil_keys.index(anchor_key)]
    stencil_bond_xyz = unit_direction(stencil_atm_xyz, stencil_anchor_xyz)

    interp_stencil_ = local_coordinate_interpreter(
        trans=true_atm_xyz,
        rot=aligning_rotation_matrix(stencil_bond_xyz, true_bond_xyz)
    )

    assert len(stencil_keys) <= len(stencil_xyzs)
    ret_xyz_dct = dict(zip(stencil_keys, map(interp_stencil_, stencil_xyzs)))

    # check that the transformation worked
    assert ret_xyz_dct[atm_key] == true_atm_xyz
    assert ret_xyz_dct[anchor_key] == true_anchor_xyz

    return ret_xyz_dct
