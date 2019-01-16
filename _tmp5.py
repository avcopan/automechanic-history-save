""" temporary script
"""
from automechanic.mol import graph3 as graph
from automechanic.mol.graph3._inchi.conn import inchi
from automechanic.mol.graph3._inchi.conn import atom_inchi_numbers
from automechanic.mol.graph3._inchi.stereo import inchi as inchi_stereo
from automechanic.mol.graph3._inchi.intco import atom_stereo_coordinates


SGR = ({0: ('C', 3, None), 1: ('C', 3, None), 2: ('C', 1, None),
        3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
        6: ('C', 2, None), 7: ('C', 1, False), 8: ('O', 0, None)},
       {frozenset({0, 2}): (1, None), frozenset({1, 3}): (1, None),
        frozenset({2, 4}): (1, False), frozenset({3, 5}): (1, False),
        frozenset({4, 6}): (1, None), frozenset({5, 7}): (1, None),
        frozenset({6, 7}): (1, None), frozenset({8, 7}): (1, None)})

SGR = graph.explicit_stereo_sites(SGR)
XYZ_DCT = atom_stereo_coordinates(SGR)

print(XYZ_DCT)
assert set(graph.atom_keys(SGR)) == set(XYZ_DCT.keys())

print(inchi(SGR))
print(atom_inchi_numbers(SGR))

SGR = graph.relabel(SGR, atom_inchi_numbers(SGR))

print(inchi(SGR))
print(atom_inchi_numbers(SGR))

print(inchi_stereo(SGR))
