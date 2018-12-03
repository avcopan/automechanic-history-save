""" pyx2z geometry interface
"""
import itertools
import numpy
import pyx2z


def resonance_count(x2m):
    """ the number of resonance structures
    """
    return x2m.resonance_count()


def bonds(x2m, ridx=0):
    """ bonds for a given resonance structure
    """
    nrncs = resonance_count(x2m)
    assert 0 <= ridx < nrncs
    natms = x2m.size()
    keys = list(map(frozenset, itertools.combinations(range(natms), r=2)))
    typs = [_bond_order(x2m, idx1, idx2) for idx1, idx2 in keys]
    bnds = {key: typ for key, typ in zip(keys, typs) if typ > 0}
    return bnds


def _bond_order(x2m, idx1, idx2, ridx=0):
    """ the bond order between two atoms
    """
    ridx = 0  # for other resonances set 0 <= ridx <= x2zms.resonance_count()
    return int(x2m.bond_order(idx1, idx2, ridx))


def from_geometry(geo):
    """ x2z molecule to geometry
    """
    x2zmg = pyx2z.MolecGeom()
    for asymb, xyz in geo:
        x2zatm = _pyx2z_atom(asymb, xyz)
        x2zmg.push_back(x2zatm)
    x2zps = pyx2z.PrimStruct(x2zmg)
    x2m = pyx2z.MolecStruct(x2zps)
    return x2m


def _pyx2z_atom(asymb, xyz):
    ang2bohr = 1.8897259886
    x2zatm = pyx2z.Atom(asymb)
    x2zatm[0], x2zatm[1], x2zatm[2] = numpy.multiply(xyz, ang2bohr)
    return x2zatm
