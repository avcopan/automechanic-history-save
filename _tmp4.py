""" temporary script
"""
from itertools import chain as _chain
from automechanic.mol import graph3 as graph
# from automechanic.mol.graph3 import coord
from automechanic.mol.graph3._dict import filter_by_value as _filter_by_value
from automechanic.mol.graph2.conn import atom_inchi_numbers


def assign_stereo_coordinates(sgr):
    """ assign stereo coordinates to this molecular graph
    """
    atm_ste_par_dct = _filter_by_value(
        graph.atom_stereo_parities(sgr), lambda val: val is not None)
    bnd_ste_par_dct = _filter_by_value(
        graph.bond_stereo_parities(sgr), lambda val: val is not None)

    atm_ste_atm_keys = list(atm_ste_par_dct.keys())
    bnd_ste_atm_keys = sorted(set(_chain(*bnd_ste_par_dct.keys())))

    ste_atm_keys = atm_ste_atm_keys + bnd_ste_atm_keys

    sgr = graph.explicit(sgr, ste_atm_keys)

    atm_ngb_keys_dct = graph.atom_neighbor_keys(sgr)
    atm_ich_num_dct = atom_inchi_numbers(sgr)

    def _assign_coordinates_recursively(atm_key,
                                        seen_atm_ngb_key=None,
                                        atm_xyz_dct=None):
        atm_ngb_keys = atm_ngb_keys_dct[atm_key]
        assert len(atm_ngb_keys) <= 4

        if atm_xyz_dct is None:
            assert seen_atm_ngb_key is None
            seen_atm_ngb_key = next(iter(atm_ngb_keys))
            atm_xyz_dct = _assign_coordinates_recursively(
                atm_key=seen_atm_ngb_key,
                seen_atm_ngb_key=atm_key,
                atm_xyz_dct={seen_atm_ngb_key: (-1, 0, 0), atm_key: (0, 0, 0)}
            )

        assert atm_key in atm_xyz_dct
        assert seen_atm_ngb_key in atm_xyz_dct
        assert seen_atm_ngb_key in atm_ngb_keys

        new_atm_ngb_keys = set(atm_ngb_keys) - {seen_atm_ngb_key}

        if atm_key in atm_ste_atm_keys:
            print("atom stereo atom key: {:d}".format(atm_key))
        elif atm_key in bnd_ste_atm_keys:
            print("bond stereo atom key: {:d}".format(atm_key))
        else:
            print('here')

        print(new_atm_ngb_keys)
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

assign_stereo_coordinates(SGR)
