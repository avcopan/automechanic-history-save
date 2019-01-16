""" temporary script
"""
from functools import partial as _partial
from itertools import chain as _chain
from automechanic.mol import graph3 as graph
from automechanic.mol.graph3._dict import filter_by_value as _filter_by_value
from automechanic.mol.graph3 import intco
from automechanic.mol.graph2.conn import atom_inchi_numbers


def stereo_coordinates(sgr):
    """ determine stereo-specific coordinates for this molecular graph
    """
    atm_par_dct = _filter_by_value(
        graph.atom_stereo_parities(sgr), lambda val: val is not None)
    bnd_par_dct = _filter_by_value(
        graph.bond_stereo_parities(sgr), lambda val: val is not None)

    atm_ste_atm_keys = list(atm_par_dct.keys())
    bnd_ste_atm_keys = sorted(set(_chain(*bnd_par_dct.keys())))
    assert sgr == graph.explicit(sgr, atm_ste_atm_keys + bnd_ste_atm_keys)

    atm_ich_num_dct = atom_inchi_numbers(sgr)
    key_sorter = _partial(sorted, key=atm_ich_num_dct.__getitem__)

    atm_ngb_keys_dct = graph.atom_neighbor_keys(sgr)

    def _recurse_coordinates(atm_key, anchor_key, xyz_dct):
        boundary_edges = ()
        if atm_key in atm_ste_atm_keys:
            atm_par = atm_par_dct[atm_key]
            xyz_dct, boundary_edges = intco.atom_stereo_coordinates(
                anchor_key, atm_key, atm_ngb_keys_dct, xyz_dct, key_sorter,
                atm_par)
        elif atm_key in bnd_ste_atm_keys:
            bnd_key, bnd_par = next(
                (bnd_key, bnd_par) for bnd_key, bnd_par in bnd_par_dct.items()
                if atm_key in bnd_key)
            xyz_dct, boundary_edges = intco.bond_stereo_coordinates(
                anchor_key, bnd_key, atm_ngb_keys_dct, xyz_dct, key_sorter,
                bnd_par)
        else:
            xyz_dct, boundary_edges = intco.nonstereo_coordinates(
                anchor_key, atm_key, atm_ngb_keys_dct, xyz_dct)

        for anchor_key_, atm_key_ in boundary_edges:
            xyz_dct = _recurse_coordinates(atm_key_, anchor_key_, xyz_dct)

        return xyz_dct

    xyz_dct = {}

    atm_keys = key_sorter(graph.atom_keys(sgr))

    if atm_keys:
        atm_key = next(iter(atm_keys))
        xyz_dct[atm_key] = (0, 0, 0)

        atm_ngb_keys = atm_ngb_keys_dct[atm_key]

        if atm_ngb_keys:
            atm_ngb_key = next(iter(atm_ngb_keys))
            xyz_dct[atm_ngb_key] = (1, 0, 0)

            xyz_dct = _recurse_coordinates(atm_ngb_key, atm_key, xyz_dct)
            xyz_dct = _recurse_coordinates(atm_key, atm_ngb_key, xyz_dct)

    return xyz_dct


SGR = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
        3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
        6: ('C', 2, None), 7: ('C', 1, False), 8: ('O', 0, None)},
       {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
        frozenset({2, 4}): (1, False), frozenset({3, 5}): (1, False),
        frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
        frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None)})

XYZ_DCT = stereo_coordinates(SGR)

print(XYZ_DCT)
print(set(graph.atom_keys(SGR)) ^ set(XYZ_DCT.keys()))

# ATM_STE_PAR_DCT = _filter_by_value(
#     graph.atom_stereo_parities(SGR), lambda val: val is not None)
# BND_STE_PAR_DCT = _filter_by_value(
#     graph.bond_stereo_parities(SGR), lambda val: val is not None)
# ATM_STE_ATM_KEYS = list(ATM_STE_PAR_DCT.keys())
# BND_STE_ATM_KEYS = sorted(set(_chain(*BND_STE_PAR_DCT.keys())))
# STE_ATM_KEYS = ATM_STE_ATM_KEYS + BND_STE_ATM_KEYS

# SGR = graph.explicit(SGR, STE_ATM_KEYS)

# ATM_ICH_NUM_DCT = atom_inchi_numbers(SGR)
# ATM_NGB_KEYS_DCT = graph.atom_neighbor_keys(SGR)

# BND_KEY = next(iter(BND_STE_PAR_DCT))
# ATM_KEY = next(iter(BND_KEY))
# ANCHOR_KEY = next(iter(ATM_NGB_KEYS_DCT[ATM_KEY]))
# KEY_SORTER = _partial(sorted, key=ATM_ICH_NUM_DCT.__getitem__)

# print(intco.bond_stereo_coordinates(
#     anchor_key=ANCHOR_KEY,
#     bnd_key=BND_KEY,
#     atm_ngb_keys_dct=ATM_NGB_KEYS_DCT,
#     xyz_dct={ANCHOR_KEY: (0, -1, 9), ATM_KEY: (0, 0, 9)},
#     key_sorter=KEY_SORTER,
#     parity=BND_STE_PAR_DCT[BND_KEY]
# ))

# print(intco.bond_stereo_coordinates(
#     anchor_key=ANCHOR_KEY,
#     bnd_key=BND_KEY,
#     atm_ngb_keys_dct=ATM_NGB_KEYS_DCT,
#     xyz_dct={ANCHOR_KEY: (0, -1, 9), ATM_KEY: (0, 0, 9)},
#     key_sorter=KEY_SORTER,
#     parity=True
# ))

# ATM_KEY = next(iter(ATM_STE_PAR_DCT))
# ANCHOR_KEY = next(iter(ATM_NGB_KEYS_DCT[ATM_KEY]))
# print(intco.atom_stereo_coordinates(
#     anchor_key=ANCHOR_KEY,
#     atm_key=ATM_KEY,
#     atm_ngb_keys_dct=ATM_NGB_KEYS_DCT,
#     xyz_dct={ANCHOR_KEY: (0, -1, 9), ATM_KEY: (0, 0, 9)},
#     key_sorter=KEY_SORTER,
#     parity=ATM_STE_PAR_DCT[ATM_KEY]
# ))

# print(intco.atom_stereo_coordinates(
#     anchor_key=ANCHOR_KEY,
#     atm_key=ATM_KEY,
#     atm_ngb_keys_dct=ATM_NGB_KEYS_DCT,
#     xyz_dct={ANCHOR_KEY: (0, -1, 9), ATM_KEY: (0, 0, 9)},
#     key_sorter=KEY_SORTER,
#     parity=True
# ))

# print(intco.nonstereo_coordinates(
#     anchor_key=ANCHOR_KEY,
#     atm_key=ATM_KEY,
#     atm_ngb_keys_dct=ATM_NGB_KEYS_DCT,
#     xyz_dct={ANCHOR_KEY: (0, -1, 9), ATM_KEY: (0, 0, 9)},
# ))
