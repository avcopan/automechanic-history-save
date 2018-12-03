""" functions operating on molecular graphs
"""
from .at import valence as _atom_valence


def atoms(gra):
    """ atomic vertices
    """
    atms, _ = gra
    return atms


def bonds(gra):
    """ bond edges
    """
    _, bnds = gra
    return bnds


def atom_count(gra):
    """ the number of atomic vertices
    """
    return len(atoms(gra))


def indices(gra):
    """ indices of the atoms in this molecular graph
    """
    return tuple(range(atom_count(gra)))


def atom_at(gra, idx):
    """ atom at this atomic vertex
    """
    assert idx in indices(gra)
    atms = atoms(gra)
    return atms[idx]


def bonds_at(gra, idx):
    """ bonds at this atomic vertex
    """
    assert idx in indices(gra)
    bnds = bonds(gra)
    atm_bnds = {key: typ for key, typ in bnds.items() if idx in key}
    return atm_bnds


def bound_electrons_at(gra, idx):
    """ bound electrons at this atomic vertex
    """
    atm_bnds = bonds_at(gra, idx)
    return sum(atm_bnds.values())


def radical_electrons_at(gra, idx):
    """ radical electrons at this atomic vertex
    """
    atm = atom_at(gra, idx)
    atm_vlc = _atom_valence(atm)
    atm_nbe = bound_electrons_at(gra, idx)
    return atm_vlc - atm_nbe


def radicals(gra):
    """ radical electrons by atom
    """
    idxs = indices(gra)
    vals = [radical_electrons_at(gra, idx) for idx in idxs]
    rads = {idx: val for idx, val in zip(idxs, vals) if val > 0}
    return rads
