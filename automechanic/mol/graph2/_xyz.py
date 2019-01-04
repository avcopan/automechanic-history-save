""" graph coordinate assignments
"""
import numpy
from ._shared import atom_keys
from ._shared import ring_keys_list
from ._shared import atom_neighbor_keys


def stereo_graph_coordinates(sgr):
    """ coordinates for this stereo graph
    """
    assert not ring_keys_list(sgr)  # currently assumes no rings
    atm_keys = atom_keys(sgr)
    atm_ngb_keys_dct = atom_neighbor_keys(sgr)
    print(atm_keys)
    print(atm_ngb_keys_dct)


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
