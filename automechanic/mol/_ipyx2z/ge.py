""" pyx2z geometry interface
"""
import numpy
import pyx2z


def bond_order(geo, idx1, idx2):
    """ the bond order between two atoms
    """
    ridx = 0  # for other resonances set 0 <= ridx <= x2zms.resonance_count()
    x2zms = _pyx2z_molec_struct(geo)
    return int(x2zms.bond_order(idx1, idx2, ridx))


def _pyx2z_molec_struct(geo):
    x2zps = _pyx2z_prim_struct(geo)
    return pyx2z.MolecStruct(x2zps)


def _pyx2z_prim_struct(geo):
    x2zmg = _pyx2z_molec_geom(geo)
    return pyx2z.PrimStruct(x2zmg)


def _pyx2z_molec_geom(geo):
    x2zmg = pyx2z.MolecGeom()
    for asymb, xyz in geo:
        x2zatm = _pyx2z_atom(asymb, xyz)
        x2zmg.push_back(x2zatm)
    return x2zmg


def _pyx2z_atom(asymb, xyz):
    ang2bohr = 1.8897259886
    x2zatm = pyx2z.Atom(asymb)
    x2zatm[0], x2zatm[1], x2zatm[2] = numpy.multiply(xyz, ang2bohr)
    return x2zatm
