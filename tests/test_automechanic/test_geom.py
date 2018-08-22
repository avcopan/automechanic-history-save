""" test the automechanic.geom module
"""
from automechanic import geom


def test__resonance_graphs():
    """ test geom.resonance_averaged_molecule_graph
    """
    mgeo = (('C', (1.04875736, -0.08263227, 0.02520069)),
            ('C', (0.35689443, -1.08220115, 0.58834534)),
            ('C', (-1.08355137, -1.09969298, 0.60845334)),
            ('C', (-1.80589237, -2.05681355, 1.20880041)),
            ('C', (-3.24001909, -1.97861000, 1.27917334)),
            ('H', (-3.71446626, -2.28637223, 2.20524241)),
            ('H', (-3.76959941, -1.26680702, 0.65403921)),
            ('H', (0.55630109, 0.75863689, -0.45190632)),
            ('H', (2.13292129, -0.08512494, 0.03758668)),
            ('H', (0.89105972, -1.90022482, 1.06515672)),
            ('H', (-1.59381524, -0.26225598, 0.13760551)),
            ('H', (-1.33092095, -2.89685079, 1.70957399)))
    mgrph1 = (('C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 7]), 1), (frozenset([9, 1]), 1),
                         (frozenset([4, 6]), 1), (frozenset([3, 4]), 1),
                         (frozenset([1, 2]), 1), (frozenset([2, 3]), 2),
                         (frozenset([3, 11]), 1), (frozenset([0, 1]), 2),
                         (frozenset([8, 0]), 1), (frozenset([10, 2]), 1),
                         (frozenset([4, 5]), 1)]))
    mgrph2 = (('C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 7]), 1), (frozenset([9, 1]), 1),
                         (frozenset([4, 6]), 1), (frozenset([1, 2]), 1),
                         (frozenset([3, 4]), 2), (frozenset([3, 11]), 1),
                         (frozenset([0, 1]), 2), (frozenset([8, 0]), 1),
                         (frozenset([2, 3]), 1), (frozenset([10, 2]), 1),
                         (frozenset([4, 5]), 1)]))
    mgrph3 = (('C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H'),
              frozenset([(frozenset([0, 7]), 1), (frozenset([1, 2]), 2),
                         (frozenset([9, 1]), 1), (frozenset([4, 6]), 1),
                         (frozenset([0, 1]), 1), (frozenset([3, 4]), 2),
                         (frozenset([3, 11]), 1), (frozenset([8, 0]), 1),
                         (frozenset([2, 3]), 1), (frozenset([10, 2]), 1),
                         (frozenset([4, 5]), 1)]))
    assert set(geom.resonance_graphs(mgeo)) == set((mgrph1, mgrph2, mgrph3))


def test__radical_sites():
    """ test geom.radical_sites
    """
    mgeo = (('C', (1.10206, 0.05263, 0.02517)),
            ('C', (2.44012, 0.03045, 0.01354)),
            ('C', (3.23570, 0.06292, 1.20436)),
            ('H', (2.86296, -0.38925, 2.11637)),
            ('H', (4.29058, 0.30031, 1.12619)),
            ('H', (0.54568, -0.01805, -0.90370)),
            ('H', (0.53167, 0.14904, 0.94292)),
            ('H', (2.97493, -0.03212, -0.93001)))
    assert geom.radical_sites(mgeo) == (0, 2)


if __name__ == '__main__':
    test__radical_sites()
    test__resonance_graphs()
