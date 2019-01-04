""" temporary script
"""
from __future__ import division
import numpy
from automechanic.mol.graph2 import stereo


def unit_integer_triple_forward(int_xyz_dir):
    """ unit integer 3-vector pointing forward along these points

    (points must be aligned along one of the axes)
    """
    assert sum(numpy.equal(int_xyz_dir, 0)) == 2

    int_xyz_forw = (int_xyz_dir / numpy.linalg.norm(int_xyz_dir)).astype(int)

    assert numpy.linalg.norm(int_xyz_forw) == 1

    return tuple(int_xyz_forw)


def unit_integer_triple_left(int_xyz_dir):
    """ unit integer 3-vector pointing left relative to this

    ('left' is an arbitrary perpendicular used the other directions)
    """
    int_xyz_forw = unit_integer_triple_forward(int_xyz_dir)

    int_x_cmp, int_y_cmp, int_z_cmp = int_xyz_forw
    if int_x_cmp != 0:
        int_xyz_left = numpy.cross(int_xyz_forw, [0, 0, -1])
    elif int_y_cmp != 0:
        int_xyz_left = numpy.cross(int_xyz_forw, [0, 0, +1])
    elif int_z_cmp != 0:
        int_xyz_left = numpy.cross(int_xyz_forw, [0, -1, 0])

    assert numpy.linalg.norm(int_xyz_left) == 1

    return tuple(int_xyz_left)


def unit_integer_triple_right(int_xyz_dir):
    """ unit integer 3-vector pointing right relative to this
    """
    int_xyz_left = unit_integer_triple_left(int_xyz_dir)
    int_xyz_right = numpy.negative(int_xyz_left)
    return tuple(int_xyz_right)


def unit_integer_triple_up(int_xyz_dir):
    """ unit integer 3-vector pointing up relative to this
    """
    int_xyz_forw = unit_integer_triple_forward(int_xyz_dir)
    int_xyz_left = unit_integer_triple_left(int_xyz_dir)
    int_xyz_up = numpy.cross(int_xyz_forw, int_xyz_left)
    return tuple(int_xyz_up)


def unit_integer_triple_down(int_xyz_dir):
    """ unit integer 3-vector pointing down relative to this
    """
    int_xyz_up = unit_integer_triple_up(int_xyz_dir)
    int_xyz_down = numpy.negative(int_xyz_up)
    return tuple(int_xyz_down)


SGR = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
        3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
        6: ('C', 2, None), 7: ('C', 1, False), 8: ('O', 0, None)},
       {frozenset({0, 2}): None, frozenset({1, 3}): None,
        frozenset({2, 4}): False, frozenset({3, 5}): False,
        frozenset({4, 6}): None, frozenset({5, 7}): None,
        frozenset({6, 7}): None, frozenset({8, 7}): None})

ATM_KEYS = stereo.atom_keys(SGR)

ATM_NGB_KEYS_DCT = stereo.atom_neighbor_keys(SGR)

ATM_PAR_DCT = stereo.atom_parities(SGR)
ATM_ICH_NUM_DCT = stereo.atom_inchi_numbers(SGR)


def _assign_neighbor_coordinates(atm_key, atm_ref_key=None):
    assert atm_key in ATM_XYZ_DCT.keys() and atm_ref_key in ATM_XYZ_DCT.keys()

    atm_ngb_keys = ATM_NGB_KEYS_DCT[atm_key]
    assert len(atm_ngb_keys) <= 4
    assert atm_ref_key in atm_ngb_keys or atm_ref_key is None

    atm_ngb_keys = sorted(set(atm_ngb_keys) - {atm_ref_key},
                          key=ATM_ICH_NUM_DCT.__getitem__)
    if atm_ngb_keys:
        atm_xyz = ATM_XYZ_DCT[atm_key]
        atm_ref_xyz = ATM_XYZ_DCT[atm_ref_key]

        int_xyz_dir = numpy.subtract(atm_xyz, atm_ref_xyz)

        if atm_key not in ATM_PAR_DCT:
            int_xyz_disps = iter((
                unit_integer_triple_forward(int_xyz_dir),
                unit_integer_triple_left(int_xyz_dir),
                unit_integer_triple_right(int_xyz_dir)))

            for atm_ngb_key, int_xyz_disp in zip(atm_ngb_keys, int_xyz_disps):
                atm_ngb_xyz = tuple(numpy.add(atm_xyz, int_xyz_disp))
                ATM_XYZ_DCT[atm_ngb_key] = atm_ngb_xyz

                _assign_neighbor_coordinates(atm_ngb_key, atm_ref_key=atm_key)
        else:
            if len(atm_ngb_keys) == 2:
                atm_par = ATM_PAR_DCT[atm_key]
                print(atm_key)
                print(atm_ref_key)
                print(atm_par)
            else:
                raise NotImplementedError


ATM_KEYS = sorted(ATM_KEYS, key=ATM_ICH_NUM_DCT.__getitem__)
ATM_KEY = ATM_KEYS[0]

ATM_NGB_KEYS = sorted(ATM_NGB_KEYS_DCT[ATM_KEY],
                      key=ATM_ICH_NUM_DCT.__getitem__)
ATM_NGB_KEY = ATM_NGB_KEYS[0]

ATM_XYZ_DCT = {}
ATM_XYZ_DCT[ATM_KEY] = (0, 0, 0)
ATM_XYZ_DCT[ATM_NGB_KEY] = (1, 0, 0)
_assign_neighbor_coordinates(ATM_NGB_KEY, atm_ref_key=ATM_KEY)
print(ATM_XYZ_DCT)
