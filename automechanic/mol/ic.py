""" functions operating on InChI strings
"""
from .ge import inchi as _inchi_from_geometry
from ._irdkit.ic import inchi as _inchi
from ._irdkit.ic import inchi_key as _inchi_key
from ._irdkit.ic import (atom_symbols_and_coordinates as
                         _atom_symbols_and_coordinates)


def inchi(ich, force_stereo=False, strict=True):
    """ recompute InChI string
    """
    _options = '-SUU' if force_stereo else ''
    ich = _inchi(ich, options=_options, strict=strict)
    return ich


def inchi_key(ich, strict=True):
    """ computes InChIKey from an InChI string
    """
    return _inchi_key(ich, strict=strict)


def geometry(ich, strict=True):
    """ cartesian geometry from an InChI string
    """
    asbs, xyzs = _atom_symbols_and_coordinates(ich, strict=strict)
    geo = tuple(zip(asbs, xyzs))
    if strict:
        ich_ = _inchi_from_geometry(geo)
        assert inchi_key(ich) == inchi_key(ich_)
    return geo
