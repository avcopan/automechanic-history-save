""" temporary script
"""
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
print(ATM_NGB_KEYS_DCT)
# XYZ_DCT = coord.assign_nonstereo_atom_neighbors(
print(coord.nonstereo_coordinates(
    atm_xyz=(0, 0, 9),
    anchor_xyz=(0, -1, 9),
    anchor_key=5,
    sorted_keys=sorted(ATM_NGB_KEYS_DCT[7], key=ATM_ICH_NUM_DCT.__getitem__)
))
print(coord.atom_stereo_coordinates(
    atm_xyz=(0, 0, 9),
    anchor_xyz=(0, -1, 9),
    anchor_key=5,
    sorted_keys=sorted(ATM_NGB_KEYS_DCT[7], key=ATM_ICH_NUM_DCT.__getitem__),
    parity=False
))
print(coord.atom_stereo_coordinates(
    atm_xyz=(0, 0, 9),
    anchor_xyz=(0, -1, 9),
    anchor_key=5,
    sorted_keys=sorted(ATM_NGB_KEYS_DCT[7], key=ATM_ICH_NUM_DCT.__getitem__),
    parity=True
))

ATM1_KEY, ATM2_KEY = next(iter(BND_STE_PAR_DCT))
ATM1_NGB_KEYS = set(ATM_NGB_KEYS_DCT[ATM1_KEY]) - {ATM2_KEY}
ATM2_NGB_KEYS = set(ATM_NGB_KEYS_DCT[ATM2_KEY]) - {ATM1_KEY}
ATM1_NGB_KEYS = sorted(ATM1_NGB_KEYS, key=ATM_ICH_NUM_DCT.__getitem__)
ATM2_NGB_KEYS = sorted(ATM2_NGB_KEYS, key=ATM_ICH_NUM_DCT.__getitem__)
ATM_KEYS = ATM1_NGB_KEYS + [ATM2_KEY] + ATM2_NGB_KEYS

print(coord.bond_stereo_coordinates(
    atm_xyz=(0, 0, 9),
    anchor_xyz=(0, -1, 9),
    anchor_key=0,
    sorted_keys=ATM_KEYS,
    parity=False
))

print(coord.bond_stereo_coordinates(
    atm_xyz=(0, 0, 9),
    anchor_xyz=(0, -1, 9),
    anchor_key=0,
    sorted_keys=ATM_KEYS,
    parity=True
))
print(ATM1_NGB_KEYS)
print(ATM2_NGB_KEYS)
