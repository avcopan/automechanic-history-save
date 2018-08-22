""" Geometry-based interface to pyx2z
"""
import numpy
import pyx2z


def graph(mgeo, res=0):
    """ molecule graph of a cartesian geometry
    """
    x2zms = _pyx2z_molec_struct(mgeo)
    natms = x2zms.size()
    atms, _ = zip(*mgeo)
    bnds = frozenset()
    for i in range(natms):
        for j in range(i):
            order = x2zms.bond_order(i, j, res)
            if order > 0:
                bnd = frozenset([i, j])
                bnds |= frozenset([(bnd, order)])
    mgrph = (atms, bnds)
    return mgrph


def number_of_resonance_graphs(mgeo):
    """ number of resonances
    """
    x2zms = _pyx2z_molec_struct(mgeo)
    nrncs = x2zms.resonance_count()
    return nrncs


def resonance_graphs(mgeo):
    """ molecule graphs of a cartesian geometry, by resonance
    """
    x2zms = _pyx2z_molec_struct(mgeo)
    natms = x2zms.size()
    nrncs = x2zms.resonance_count()
    atms, _ = zip(*mgeo)
    mgrphs = []
    for ridx in range(nrncs):
        bnds = frozenset()
        for i in range(natms):
            for j in range(i):
                order = x2zms.bond_order(i, j, ridx)
                if order > 0:
                    bnd = frozenset([i, j])
                    bnds |= frozenset([(bnd, order)])
        mgrphs.append((atms, bnds))
    return tuple(mgrphs)


def radical_sites(mgeo):
    """ radical sites of a molecule
    """
    x2zms = _pyx2z_molec_struct(mgeo)
    natms = x2zms.size()
    idxs = tuple(i for i in range(natms) if x2zms.is_radical(i))
    return idxs


def _pyx2z_molec_struct(mgeo):
    x2zps = _pyx2z_prim_struct(mgeo)
    return pyx2z.MolecStruct(x2zps)


def _pyx2z_prim_struct(mgeo):
    x2zmg = _pyx2z_molec_geom(mgeo)
    return pyx2z.PrimStruct(x2zmg)


def _pyx2z_molec_geom(mgeo):
    x2zmg = pyx2z.MolecGeom()
    for asymb, xyz in mgeo:
        x2zatm = _pyx2z_atom(asymb, xyz)
        x2zmg.push_back(x2zatm)
    return x2zmg


def _pyx2z_atom(asymb, xyz):
    ang2bohr = 1.8897259886
    x2zatm = pyx2z.Atom(asymb)
    x2zatm[0], x2zatm[1], x2zatm[2] = numpy.multiply(xyz, ang2bohr)
    return x2zatm
