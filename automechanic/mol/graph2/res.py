""" resonance graph library

vertices: atomic symbols and implicit hydrogen counts
    (('O', 1), ('C', 1), ('C', 1), ...)
edges: bond connectivity and bond orders
    {{0, 1}: 1, {1, 2}: 2, ...}
"""
import numpy
from .base import vertex_keys as _vertex_keys
from .base import vertex_edges as _vertex_edges
from .base import edge_keys as _edge_keys
from ._shared import atomic_symbols
from ._shared import hydrogen_counts
from ._shared import valences
from ._shared import make_hydrogens_implicit
from ..atom import lone_pair_count as _atom_lone_pair_count
from .._molfile import from_data as _mlf_from_data
from .._molfile import to_inchi as _mlf_to_inchi
from .._inchi_aux import numbering as _ich_aux_numbering


def bond_counts(rgr):
    """ total bond counts for each atom, accounting for bond orders
    """
    hcnts = hydrogen_counts(rgr)
    xcnts = tuple(sum(_vertex_edges(rgr, key).values())
                  for key in _vertex_keys(rgr))
    bnd_cnts = numpy.add(hcnts, xcnts)
    return bnd_cnts


def radical_electron_counts(rgr):
    """ the number of radical electrons for each atom
    """
    rad_elec_cnts = numpy.subtract(valences(rgr), bond_counts(rgr))
    # assumes no hypervalent atoms, so check
    assert all(numpy.greater_equal(rad_elec_cnts, 0))
    return tuple(rad_elec_cnts)


def radical_atom_keys(rgr):
    """ keys to atoms with radical electrons
    """
    return tuple(key for key, cnt in enumerate(radical_electron_counts(rgr))
                 if cnt > 0)


def open_pi_bond_keys(rgr):
    """ open pi bond (adjacent radical site) keys
    """
    bnd_keys = _edge_keys(rgr)
    rad_keys = set(radical_atom_keys(rgr))
    return tuple(bnd_key for bnd_key in bnd_keys if bnd_key <= rad_keys)


def lone_pair_counts(cgr):
    """ lone pair counts for each atom (independent of connectivity)
    """
    syms = atomic_symbols(cgr)
    lp_cnts = tuple(map(_atom_lone_pair_count, syms))
    return lp_cnts


def sigma_bond_counts(rgr):
    """ sigma bond counts for each atm
    """
    hcnts = hydrogen_counts(rgr)
    xcnts = tuple(len(_vertex_edges(rgr, key)) for key in _vertex_keys(rgr))
    sig_bnd_cnts = numpy.add(hcnts, xcnts)
    return sig_bnd_cnts


def pi_bond_forming_resonances(rgr):
    """ recursively determine all pi-bond forming resonances
    """
    print(rgr)


def possible_spin_multiplicities(rgr):
    """ possible spin multiplicities for this resonance graph
    """
    rad_elec_cnt = sum(radical_electron_counts(rgr))
    return _spin_multiplicities(rad_elec_cnt)


def _spin_multiplicities(nelecs):
    assert nelecs >= 0
    mult_max = nelecs + 1
    mult_min = 2 if (mult_max % 2 == 0) else 1
    mults = tuple(range(mult_min, mult_max+1, 2))
    return mults


def inchi(rgr):
    """ InChI string of a connectivity graph
    """
    ich, _ = _inchi_with_numbering(rgr)
    return ich


def inchi_numbering(rgr):
    """ InChI numbering of backbone atoms
    """
    _, nums = _inchi_with_numbering(rgr)
    return nums


def _inchi_with_numbering(rgr):
    ich_nums_hardcoded = _inchi_with_numbering_hardcoded_exception(rgr)
    if ich_nums_hardcoded is not None:
        # MOLFile strings can't handle atoms with more than 2 radical
        # electrons, so we hardcode some special cases like the carbon atom.
        ich, nums = ich_nums_hardcoded
    else:
        # get the InChI string from its MOLFile string
        bnd_keys = _edge_keys(rgr)
        atm_syms = atomic_symbols(rgr)
        atm_rad_cnts = radical_electron_counts(rgr)
        atm_bnd_cnts = bond_counts(rgr)

        mlf = _mlf_from_data(atm_syms=atm_syms,
                             atm_bnd_cnts=atm_bnd_cnts,
                             atm_rad_cnts=atm_rad_cnts,
                             bnd_keys=bnd_keys)
        ich, ich_aux = _mlf_to_inchi(mlf, with_aux_info=True)
        nums = _ich_aux_numbering(ich_aux)
    return ich, nums


def _inchi_with_numbering_hardcoded_exception(rgr):
    ret = None
    if rgr == ((('C', 0),), {}):
        ret = ('InChI=1S/C', (0,))
    elif rgr == ((('N', 0),), {}):
        ret = ('InChI=1S/N', (0,))
    elif make_hydrogens_implicit(rgr) == ((('C', 1),), {}):
        nums = (1,) if atomic_symbols(rgr) == ('H', 'C') else (0,)
        ret = ('InChI=1S/CH/h1H', nums)
    return ret
