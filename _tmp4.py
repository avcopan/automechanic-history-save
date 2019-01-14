""" temporary script
"""
import numpy
from functools import partial as _partial
# from itertools import chain as _chain
from automechanic.mol import graph3 as graph
from automechanic.mol.graph3 import coord
from automechanic.mol.graph3._dict import filter_by_value as _filter_by_value
from automechanic.mol.graph2.conn import atom_inchi_numbers


def assign_stereo_coordinates(sgr):
    """ assign stereo coordinates to this molecular graph
    """
    atm_ste_par_dct = _filter_by_value(
        graph.atom_stereo_parities(sgr), lambda val: val is not None)
    # bnd_ste_par_dct = _filter_by_value(
    #     graph.bond_stereo_parities(sgr), lambda val: val is not None)

    atm_ste_atm_keys = list(atm_ste_par_dct.keys())
    bnd_ste_atm_keys = []
    # bnd_ste_atm_keys = sorted(set(_chain(*bnd_ste_par_dct.keys())))

    ste_atm_keys = atm_ste_atm_keys + bnd_ste_atm_keys

    sgr = graph.explicit(sgr, ste_atm_keys)

    atm_ngb_keys_dct = graph.atom_neighbor_keys(sgr)
    atm_ich_num_dct = atom_inchi_numbers(sgr)

    def _assign_coordinates_recursively(atm_key,
                                        anchor_atm_ngb_key=None,
                                        atm_xyz_dct=None):
        atm_ngb_keys = atm_ngb_keys_dct[atm_key]
        assert len(atm_ngb_keys) <= 4

        if atm_xyz_dct is None:
            assert anchor_atm_ngb_key is None
            anchor_atm_ngb_key = next(iter(atm_ngb_keys))
            atm_xyz_dct = _assign_coordinates_recursively(
                atm_key=anchor_atm_ngb_key,
                anchor_atm_ngb_key=atm_key,
                atm_xyz_dct={anchor_atm_ngb_key: (-1, 0, 0),
                             atm_key: (0, 0, 0)}
            )

        assert atm_key in atm_xyz_dct
        assert anchor_atm_ngb_key in atm_xyz_dct
        assert anchor_atm_ngb_key in atm_ngb_keys

        atm_ngb_keys = sorted(atm_ngb_keys, key=atm_ich_num_dct.__getitem__)

        print(atm_xyz_dct)
        atm_xyz = atm_xyz_dct[atm_key]
        anchor_pos = atm_ngb_keys.index(anchor_atm_ngb_key)
        anchor_xyz = numpy.subtract(atm_xyz_dct[anchor_atm_ngb_key], atm_xyz)

        if atm_key in atm_ste_atm_keys:
            print("atom stereo atom key: {:d}".format(atm_key))
        elif atm_key in bnd_ste_atm_keys:
            print("bond stereo atom key: {:d}".format(atm_key))
        else:
            local_xyzs = coord.lib.NON_STEREO_ATOM_NEIGHBORS
            local_anchor_xyz = local_xyzs[anchor_pos]
            align_ = coord.rot.aligning_rotator(local_anchor_xyz, anchor_xyz)
            trans_ = _partial(numpy.add, atm_xyz)
            xyzs = list(map(tuple, map(trans_, map(align_, local_xyzs))))

            atm_ngb_keys.pop(anchor_pos)
            xyzs.pop(anchor_pos)
            assert len(xyzs) > len(atm_ngb_keys)
            atm_xyz_dct.update(dict(zip(atm_ngb_keys, xyzs)))

            for atm_ngb_key in atm_ngb_keys:
                atm_xyz_dct = _assign_coordinates_recursively(
                    atm_key=atm_ngb_key,
                    anchor_atm_ngb_key=atm_key,
                    atm_xyz_dct=atm_xyz_dct
                )

        return atm_xyz_dct

    start_atm_key = sorted(graph.atom_keys(sgr),
                           key=atm_ich_num_dct.__getitem__)[0]
    _assign_coordinates_recursively(start_atm_key)


SGR = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
        3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
        6: ('C', 2, None), 7: ('C', 1, False), 8: ('O', 0, None)},
       {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
        frozenset({2, 4}): (1, False), frozenset({3, 5}): (1, False),
        frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
        frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None)})

atm_xyz_dct = assign_stereo_coordinates(SGR)
print(atm_xyz_dct)
