""" functions operating on cartesian geometries (angstroms)
"""
import numpy
from ._ipyx2z import from_geometry as _x2m_from_geometry
from ._ipyx2z import bonds as _x2m_bonds
from ._irdkit import from_molfile as _rdm_from_molfile2
from ._irdkit import to_inchi as _rdm_to_inchi
from .graph import edges as _graph_edges
from .graph.res import radical_sites as _resonance_graph_radical_sites
from ..rere.pattern import escape as _escape
from ..rere.pattern import capturing as _capturing
from ..rere.pattern import zero_or_more as _zero_or_more
from ..rere.pattern_lib import UNSIGNED_INTEGER as _UNSIGNED_INTEGER
from ..rere.find import split as _split
from ..rere.find import single_capture as _single_capture


def inchi(geo):
    """ InChI string of a cartesian geometry
    """
    ich, _ = inchi_with_order(geo)
    return ich


def inchi_with_order(geo):
    """ InChI string of a cartesian geometry
    """
    mbl = molfile2(geo)
    rdm = _rdm_from_molfile2(mbl)
    ich, ich_aux = _rdm_to_inchi(rdm, with_aux_info=True)
    ich_ord = _parse_inchi_order_from_auxinfo(ich_aux)
    return ich, ich_ord


def _parse_inchi_order_from_auxinfo(ich_aux):
    _comma = _escape(',')
    _pattern = _escape('/N:') + _capturing(
        _zero_or_more(_UNSIGNED_INTEGER + _comma) + _UNSIGNED_INTEGER)
    one_index_order_str = _single_capture(_pattern, ich_aux)
    one_index_order = tuple(map(int, _split(_comma, one_index_order_str)))
    order = tuple(numpy.subtract(one_index_order, 1))
    return order


def atoms(geo):
    """ molecule atomic symbols
    """
    atms, _ = zip(*geo)
    return atms


def resonance_graph(geo):
    """ molecule resonance graph
    """
    asbs, _ = zip(*geo)
    bnds = {}
    if len(geo) > 1:
        x2m = _x2m_from_geometry(geo)
        bnds = _x2m_bonds(x2m, ridx=0)

    return (asbs, bnds)


def molfile2(geo):
    """ mol block string of a cartesian geometry
    """
    sections = []
    # form the header:
    _head = ('\nautomechanic\n\n'
             '{natms:>3d}{nbnds:>3d}  0  0  0  0  0  0  0  0999 V2000')
    rgr = resonance_graph(geo)
    bnds = _graph_edges(rgr)
    natms = len(geo)
    nbnds = len(bnds)
    head = _head.format(natms=natms, nbnds=nbnds)
    sections.append(head)

    # form the atoms block of the body:
    _atm_line = ('{x:>10.4f}{y:>10.4f}{z:>10.4f} {a:<3s}'
                 ' 0  0  0  0  0  0  0  0  0  0  0  0')
    asbs, xyzs = zip(*geo)
    body_atms = '\n'.join(
        _atm_line.format(x=x, y=y, z=z, a=a)
        for a, (x, y, z) in zip(asbs, xyzs))
    sections.append(body_atms)

    # form the bonds block of the body:
    if bnds:
        _bnd_line = '{i:>3d}{j:>3d}{t:>3d}  0  0  0  0'
        bkeys, btyps = zip(*bnds.items())
        bkeys = list(map(sorted, bkeys))
        one_index_bkeys = numpy.add(bkeys, 1)
        body_bnds = '\n'.join(
            _bnd_line.format(i=i, j=j, t=t)
            for (i, j), t in zip(one_index_bkeys, btyps))
        sections.append(body_bnds)

    # form the properties block of the body:
    _rad_line = 'M  RAD  1{i:>4d}{m:>4d}'
    rads = _resonance_graph_radical_sites(rgr)
    if rads:
        rkeys, rvals = zip(*rads.items())
        mults = numpy.add(rvals, 1)
        one_index_rkeys = numpy.add(rkeys, 1)
        body_rads = '\n'.join(
            _rad_line.format(i=i, m=m)
            for i, m in zip(one_index_rkeys, mults))
        sections.append(body_rads)

    # add the footer
    foot = 'M  END'
    sections.append(foot)

    mbl = '\n'.join(sections)
    return mbl
