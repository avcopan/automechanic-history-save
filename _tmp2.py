""" temporary script
"""
import numpy
from automechanic.mol.graph2 import stereo


def perpendicular_unit_integer_triple(int_xyz_unit):
    """ perpendicular to an integer unit 3-vector
    """
    int_xyz_unit = numpy.array(int_xyz_unit)

    assert issubclass(int_xyz_unit.dtype.type, numpy.integer)
    assert numpy.linalg.norm(int_xyz_unit) == 1

    int_x_cmp, int_y_cmp, int_z_cmp = int_xyz_unit
    if int_x_cmp != 0:
        int_xyz_perp = numpy.cross(int_xyz_unit, [0, 0, -1])
    elif int_y_cmp != 0:
        int_xyz_perp = numpy.cross(int_xyz_unit, [0, 0, +1])
    elif int_z_cmp != 0:
        int_xyz_perp = numpy.cross(int_xyz_unit, [0, -1, 0])

    assert numpy.linalg.norm(int_xyz_perp) == 1
    return tuple(int_xyz_perp)


SGR = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
        3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
        6: ('C', 2, None), 7: ('C', 1, False), 8: ('O', 0, None)},
       {frozenset({0, 2}): None, frozenset({1, 3}): None,
        frozenset({2, 4}): False, frozenset({3, 5}): False,
        frozenset({4, 6}): None, frozenset({5, 7}): None,
        frozenset({6, 7}): None, frozenset({8, 7}): None})

ATM_KEYS = stereo.atom_keys(SGR)

ATM_XYZ_DCT = {}
ATM_XYZ_DCT[0] = (0, 0, 0)
ATM_XYZ_DCT[2] = (1, 0, 0)

ATM_NGB_KEYS_DCT = stereo.atom_neighbor_keys(SGR)


def _assign_neighbor_coordinates(atm_key, atm_ref_key):
    assert atm_key in ATM_XYZ_DCT.keys() and atm_ref_key in ATM_XYZ_DCT.keys()

    atm_ngb_keys = ATM_NGB_KEYS_DCT[atm_key]
    assert len(atm_ngb_keys) <= 4
    assert atm_ref_key in atm_ngb_keys

    atm_ngb_keys = sorted(set(atm_ngb_keys) - {atm_ref_key})
    if atm_ngb_keys:
        atm_xyz = ATM_XYZ_DCT[atm_key]
        atm_ref_xyz = ATM_XYZ_DCT[atm_ref_key]

        disp1_xyz = tuple(numpy.subtract(atm_xyz, atm_ref_xyz))
        disp2_xyz = perpendicular_unit_integer_triple(disp1_xyz)
        disp3_xyz = tuple(numpy.negative(disp2_xyz))

        disp_xyzs = iter((disp1_xyz, disp2_xyz, disp3_xyz))
        for atm_ngb_key, disp_xyz in zip(atm_ngb_keys, disp_xyzs):
            atm_ngb_xyz = tuple(numpy.add(atm_xyz, disp_xyz))
            ATM_XYZ_DCT[atm_ngb_key] = atm_ngb_xyz

            _assign_neighbor_coordinates(atm_ngb_key, atm_ref_key=atm_key)


_assign_neighbor_coordinates(2, atm_ref_key=0)
print(ATM_XYZ_DCT)
