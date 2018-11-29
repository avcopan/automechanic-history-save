""" functions operating on cartesian geometries (angstroms)
"""
import itertools
import numpy
from ._ipyx2z.ge import bond_order as _bond_order
from ._irdkit.mb import (inchi_with_auxinfo as
                         _inchi_with_auxinfo_from_mol_block)
from ..rere.pattern import escape as _escape
from ..rere.pattern import capturing as _capturing
from ..rere.pattern import zero_or_more as _zero_or_more
from ..rere.pattern_lib import UNSIGNED_INTEGER as _UNSIGNED_INTEGER
from ..rere.find import split as _split
from ..rere.find import single_capture as _single_capture


def inchi(geo, force_stereo=False, strict=True):
    """ InChI string of a cartesian geometry
    """
    ich, _ = inchi_with_order(geo, force_stereo=force_stereo, strict=strict)
    return ich


def inchi_with_order(geo, force_stereo=False, strict=True):
    """ InChI string of a cartesian geometry
    """
    _mbl = mol_block(geo)
    _options = '-SUU' if force_stereo else ''
    ich, ich_aux = _inchi_with_auxinfo_from_mol_block(
        _mbl, options=_options, strict=strict)
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


def bonds(geo):
    """ molecule bonds, as a dictionary
    """
    natms = len(geo)
    keys = list(map(frozenset, itertools.combinations(range(natms), r=2)))
    typs = [_bond_order(geo, idx1, idx2) for idx1, idx2 in keys]
    bnds = {key: typ for key, typ in zip(keys, typs) if typ > 0}
    return bnds


def mol_block(geo):
    """ mol block string of a cartesian geometry
    """
    _head = ('\nautomechanic\n\n'
             '{natms:>3d}{nbnds:>3d}  0  0  0  0  0  0  0  0999 V2000')
    _atm_line = ('{x:>10.4f}{y:>10.4f}{z:>10.4f} {a:<3s}'
                 ' 0  0  0  0  0  0  0  0  0  0  0  0')
    _bnd_line = '{i:>3d}{j:>3d}{t:>3d}  0  0  0  0'
    foot = 'M  END'

    # form the atoms block of the body:
    asbs, xyzs = zip(*geo)
    body_atms = '\n'.join(
        _atm_line.format(x=x, y=y, z=z, a=a)
        for a, (x, y, z) in zip(asbs, xyzs))

    # form the header:
    bnds = bonds(geo)
    natms = len(geo)
    nbnds = len(bnds)
    head = _head.format(natms=natms, nbnds=nbnds)

    # form the bonds block of the body:
    bkeys, btyps = zip(*bnds.items())
    bkeys = list(map(sorted, bkeys))
    one_index_bkeys = numpy.add(bkeys, 1)
    body_bnds = '\n'.join(
        _bnd_line.format(i=i, j=j, t=t)
        for (i, j), t in zip(one_index_bkeys, btyps))

    mbl = '\n'.join([head, body_atms, body_bnds, foot])
    return mbl
