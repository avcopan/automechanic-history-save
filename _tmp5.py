""" temporary script
"""
from functools import partial as _partial
from itertools import chain as _chain
from automechanic.mol import graph3 as graph
from automechanic.mol.graph3._dict import filter_by_value as _filter_by_value
from automechanic.mol.graph3 import coord
from automechanic.mol.graph2.conn import atom_inchi_numbers


SGR = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
        3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
        6: ('C', 2, None), 7: ('C', 1, False), 8: ('O', 0, None)},
       {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
        frozenset({2, 4}): (1, False), frozenset({3, 5}): (1, False),
        frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
        frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None)})

ATM_STE_PAR_DCT = _filter_by_value(
    graph.atom_stereo_parities(SGR), lambda val: val is not None)
BND_STE_PAR_DCT = _filter_by_value(
    graph.bond_stereo_parities(SGR), lambda val: val is not None)
ATM_STE_ATM_KEYS = list(ATM_STE_PAR_DCT.keys())
BND_STE_ATM_KEYS = sorted(set(_chain(*BND_STE_PAR_DCT.keys())))
STE_ATM_KEYS = ATM_STE_ATM_KEYS + BND_STE_ATM_KEYS

print(BND_STE_ATM_KEYS)

SGR = graph.explicit(SGR, STE_ATM_KEYS)

ATM_ICH_NUM_DCT = atom_inchi_numbers(SGR)
ATM_NGB_KEYS_DCT = graph.atom_neighbor_keys(SGR)
# XYZ_DCT = coord.assign_nonstereo_atom_neighbors(
print(coord.nonstereo_coordinates(
    atm_xyz=(0, 0, 9),
    anchor_xyz=(0, -1, 9),
    anchor_key=5,
    sorted_keys=sorted(ATM_NGB_KEYS_DCT[7], key=ATM_ICH_NUM_DCT.__getitem__)
))
# print(coord.atom_stereo_coordinates(
#     atm_xyz=(0, 0, 9),
#     anchor_xyz=(0, -1, 9),
#     anchor_key=5,
#     sorted_keys=sorted(ATM_NGB_KEYS_DCT[7], key=ATM_ICH_NUM_DCT.__getitem__),
#     parity=False
# ))
# print(coord.atom_stereo_coordinates(
#     atm_xyz=(0, 0, 9),
#     anchor_xyz=(0, -1, 9),
#     anchor_key=5,
#     sorted_keys=sorted(ATM_NGB_KEYS_DCT[7], key=ATM_ICH_NUM_DCT.__getitem__),
#     parity=True
# ))

BND_KEY = next(iter(BND_STE_PAR_DCT))
ATM_KEY = next(iter(BND_KEY))
ANCHOR_KEY = next(iter(ATM_NGB_KEYS_DCT[ATM_KEY]))
KEY_SORTER = _partial(sorted, key=ATM_ICH_NUM_DCT.__getitem__)

print(coord.bond_stereo_coordinates(
    anchor_key=ANCHOR_KEY,
    bnd_key=BND_KEY,
    atm_ngb_keys_dct=ATM_NGB_KEYS_DCT,
    xyz_dct={ANCHOR_KEY: (0, -1, 9), ATM_KEY: (0, 0, 9)},
    key_sorter=KEY_SORTER,
    parity=BND_STE_PAR_DCT[BND_KEY]
))

print(coord.bond_stereo_coordinates(
    anchor_key=ANCHOR_KEY,
    bnd_key=BND_KEY,
    atm_ngb_keys_dct=ATM_NGB_KEYS_DCT,
    xyz_dct={ANCHOR_KEY: (0, -1, 9), ATM_KEY: (0, 0, 9)},
    key_sorter=KEY_SORTER,
    parity=True
))

ATM_KEY = next(iter(ATM_STE_PAR_DCT))
ANCHOR_KEY = next(iter(ATM_NGB_KEYS_DCT[ATM_KEY]))
print(coord.atom_stereo_coordinates(
    anchor_key=ANCHOR_KEY,
    atm_key=ATM_KEY,
    atm_ngb_keys_dct=ATM_NGB_KEYS_DCT,
    xyz_dct={ANCHOR_KEY: (0, -1, 9), ATM_KEY: (0, 0, 9)},
    key_sorter=KEY_SORTER,
    parity=ATM_STE_PAR_DCT[ATM_KEY]
))
