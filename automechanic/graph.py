""" graph functions
"""
from functools import partial
from itertools import permutations
from more_itertools import unique_everseen
import numpy
from .atom import valence
from .perm import compose
from .perm import inverse
from .timeout import timeout


def atoms(mgrph):
    """ the atoms in this molecular graph
    """
    atms, _ = mgrph
    return atms


def formula(mgrph):
    """ the molecular formula of this molecular graph
    """
    atms = atoms(mgrph)
    return dict(_formula(atms))


def number_of_heavy_atoms(mgrph):
    """ the number of heavy atoms in this molecular graph
    """
    fml = formula(mgrph)
    return sum(fml[atm] for atm in fml if atm != 'H')


def bonds(mgrph):
    """ the bonds in this molecular graph
    """
    _, bnds = mgrph
    return bnds


def bond_keys(mgrph):
    """ the bond keys in this molecular graph
    """
    bnds = bonds(mgrph)
    bkeys, _ = _bond_keys_and_types(bnds)
    return bkeys


def bond_orders(mgrph):
    """ bond orders returned as a dictionary
    """
    bnds = bonds(mgrph)
    bkeys, btyps = _bond_keys_and_types(bnds)
    bdct = dict(zip(bkeys, btyps))
    return bdct


def radical_sites(mgrph):
    """ radical sites of a molecular graph
    """
    natms = len(atoms(mgrph))
    idxs = tuple(idx for idx in range(natms)
                 if atom_free_electrons(mgrph, idx) > 0)
    return idxs


def neighborhood_formulas(mgrph, order=1):
    """ map atoms to a list of neighborhood formulas
    """
    natms = len(atoms(mgrph))
    atm_nbhd_fml_ = partial(atom_neighborhood_formula, mgrph, order=order)
    nbhd_fmls = tuple(map(atm_nbhd_fml_, range(natms)))
    return nbhd_fmls


def atom_bonds(mgrph, idx):
    """ the bonds of an atom in the molecule graph
    """
    bnds = bonds(mgrph)
    atm_bnds = frozenset([(bkey, btyp) for bkey, btyp in bnds if idx in bkey])
    return atm_bnds


def delete_atom(mgrph, idx):
    """ delete an atom from the molecule graph
    """
    atms = atoms(mgrph)
    bnds = bonds(mgrph)
    natms = len(atms)

    keep_idxs = [i for i in range(natms) if i != idx]
    ret_atms = tuple(atms[i] for i in keep_idxs)

    atm_bnds = atom_bonds(mgrph, idx)
    keep_bnds = bnds - atm_bnds
    idx_dct = dict(map(reversed, enumerate(keep_idxs)))
    ret_bnds = _change_indices(keep_bnds, idx_dct)
    ret_mgrph = (ret_atms, ret_bnds)
    return ret_mgrph


def bind_atom(mgrph, idx, atm, order=1):
    """ bind an atom at a given site in the molecule graph
    """
    atms = atoms(mgrph)
    bnds = bonds(mgrph)
    ret_atms = atms + (atm,)
    other_idx = len(atms)
    bnd = (frozenset([idx, other_idx]), order)
    ret_bnds = set(bnds)
    ret_bnds.add(bnd)
    ret_bnds = frozenset(ret_bnds)
    ret_mgrph = (ret_atms, ret_bnds)
    return ret_mgrph


def atom_neighborhood_indices(mgrph, idx):
    """ atoms adjacent to this one in a molecule graph, by index
    """
    atm_bnds = atom_bonds(mgrph, idx)
    atm_bkeys, _ = _bond_keys_and_types(atm_bnds)
    atm_scnds = tuple(_other_vertex(bkey, idx) for bkey in atm_bkeys)
    return frozenset(atm_scnds)


def atom_neighborhood_formula(mgrph, idx, order=1):
    """ atoms adjacent to this one in a molecule graph
    """
    assert order == int(order) and order >= 1
    order = int(order)
    atms = atoms(mgrph)
    nbhd_idxs = atom_neighborhood_indices(mgrph, idx)
    if order == 1:
        nbhd_atms = tuple(map(atms.__getitem__, nbhd_idxs))
        nbhd_fml = _formula(nbhd_atms)
    elif order > 1:
        nbhd_fml_ = partial(atom_neighborhood_formula, mgrph, order=order-1)
        nbhd_items = tuple(map(nbhd_fml_, nbhd_idxs))
        nbhd_fml = _formula(nbhd_items)
    return nbhd_fml


def atom_bond_count(mgrph, idx):
    """ bond count of an atom in a molecule graph
    """
    atm_bnds = atom_bonds(mgrph, idx)
    _, atm_btyps = _bond_keys_and_types(atm_bnds)
    return sum(atm_btyps)


def atom_free_electrons(mgrph, idx):
    """ number of unbound valence electrons for an atom in a molecular graph
    """
    atms = atoms(mgrph)
    vlnc = valence(atms[idx])
    bcnt = atom_bond_count(mgrph, idx)
    return vlnc - bcnt


def disjoint_union(mgrph1, mgrph2):
    """ combine two graphs into one
    """
    atms = atoms(mgrph1) + atoms(mgrph2)
    natms1 = len(atoms(mgrph1))
    natms2 = len(atoms(mgrph2))
    bnds2 = bonds(mgrph2)
    bnds2_idx_dct = dict(enumerate(range(natms1, natms1 + natms2)))
    bnds = bonds(mgrph1) | _change_indices(bnds2, bnds2_idx_dct)
    mgrph = (atms, bnds)
    return mgrph


def increment_bond_order(mgrph, idx1, idx2, incr=+1):
    """ increment the order of a bond
    """
    assert incr <= atom_free_electrons(mgrph, idx1)
    assert incr <= atom_free_electrons(mgrph, idx2)
    atms = atoms(mgrph)
    bdct = bond_orders(mgrph)
    bkey = frozenset([idx1, idx2])
    btyp = bdct[bkey] if bkey in bdct else 0
    assert btyp + incr >= 0
    bdct[bkey] = btyp + incr
    bnds = frozenset(bdct.items())
    mgrph = (atms, bnds)
    return mgrph


def bind(mgrph1, mgrph2, idx1, idx2, order=1):
    """ bind molecule graphs
    """
    shifted_idx2 = idx2 + len(atoms(mgrph1))
    mgrph = disjoint_union(mgrph1, mgrph2)
    mgrph = increment_bond_order(mgrph, idx1, shifted_idx2, incr=order)
    return mgrph


def _other_vertex(bkey, idx):
    ret_idx, = bkey - frozenset([idx])
    return ret_idx


def _shift_bond_keys(bnds, shift):
    if not bnds:
        ret_bnds = frozenset()
    else:
        bkeys, btyps = _bond_keys_and_types(bnds)
        frsts, scnds = _unzip_bond_keys(bkeys)
        ret_frsts = numpy.add(frsts, shift)
        ret_scnds = numpy.add(scnds, shift)
        ret_bkeys = map(frozenset, zip(ret_frsts, ret_scnds))
        ret_bnds = frozenset(zip(ret_bkeys, btyps))
    return ret_bnds


def _bond_keys_and_types(bnds):
    bkeys, btyps = zip(*bnds) if bnds else ((), ())
    return bkeys, btyps


def _unzip_bond_keys(bkeys):
    frsts, scnds = zip(*bkeys) if bkeys else ((), ())
    return frsts, scnds


def isomorphic(mgrph1, mgrph2):
    """ are these graphs isomorphic?
    """
    return bool(isomorphism(mgrph1, mgrph2))


def isomorphism(mgrph, target_mgrph):
    """ isomorphism between graphs, if there is one
    """
    return _isomorphism3(mgrph, target_mgrph)


def sort_by_atoms(mgrph):
    """ sort the atoms in a molecular graph by atomic symbol
    """
    atms, _ = mgrph
    pmt = tuple(numpy.argsort(atms))
    ret_mgrph = permute_atoms(mgrph, pmt)
    return ret_mgrph


def permute_atoms(mgrph, pmt):
    """ permute the atoms in a molecular graph
    """
    atms, bnds = mgrph
    ret_atms = tuple(atms[i] for i in pmt)
    inv_pmt = inverse(pmt)
    idx_dct = dict(enumerate(inv_pmt))
    ret_bnds = _change_indices(bnds, idx_dct)
    return ret_atms, ret_bnds


@timeout(20)
def _isomorphism1(mgrph, target_mgrph):
    """ are these graphs isomorphic? (slow)
    """
    iso = None
    if formula(mgrph) == formula(target_mgrph):
        natms = len(atoms(target_mgrph))
        for pmt in permutations(range(natms)):
            pmt_mgrph = permute_atoms(mgrph, pmt)
            if pmt_mgrph == target_mgrph:
                iso = pmt
                break
    return iso


@timeout(20)
def _isomorphism2(mgrph, target_mgrph):
    """ are these graphs isomorphic? (slightly faster)
    """
    iso = None
    if formula(mgrph) == formula(target_mgrph):
        r_srt = numpy.argsort(atoms(mgrph))
        t_srt = numpy.argsort(atoms(target_mgrph))
        r_mgrph = permute_atoms(mgrph, r_srt)
        t_mgrph = permute_atoms(target_mgrph, t_srt)

        atms = atoms(r_mgrph)
        idxs_by_atm = [[i for i, a in enumerate(atms) if a == atm]
                       for atm in unique_everseen(atms)]

        for pmt in _flat_product_permutations(*idxs_by_atm):
            p_mgrph = permute_atoms(r_mgrph, pmt)
            if p_mgrph == t_mgrph:
                srt_iso = pmt
                iso = compose(inverse(t_srt), srt_iso, r_srt)
                break
    return iso


def _argsort(seq):
    return tuple(sorted(range(len(seq)), key=seq.__getitem__))


@timeout(20)
def _isomorphism3(mgrph, target_mgrph, max_order=15):
    """ are these graphs isomorphic? (slightly faster)
    """
    iso = None
    if formula(mgrph) == formula(target_mgrph):
        nhvys = number_of_heavy_atoms(mgrph)
        order = min(nhvys + 1, max_order)
        r_nbhd_fmls = neighborhood_formulas(mgrph, order=order)
        t_nbhd_fmls = neighborhood_formulas(target_mgrph, order=order)
        if sorted(r_nbhd_fmls) == sorted(t_nbhd_fmls):
            r_srt = _argsort(r_nbhd_fmls)
            t_srt = _argsort(t_nbhd_fmls)
            r_mgrph = permute_atoms(mgrph, r_srt)
            t_mgrph = permute_atoms(target_mgrph, t_srt)

            nbhd_fmls = tuple(r_nbhd_fmls[i] for i in r_srt)
            idxs_by_nbhd = [[i for i, f in enumerate(nbhd_fmls) if f == fml]
                            for fml in unique_everseen(nbhd_fmls)]

            for pmt in _flat_product_permutations(*idxs_by_nbhd):
                p_mgrph = permute_atoms(r_mgrph, pmt)
                if p_mgrph == t_mgrph:
                    srt_iso = pmt
                    iso = compose(inverse(t_srt), srt_iso, r_srt)
                    break
    return iso


def _flat_product_permutations(*seqs):
    perm_creators = tuple(map(_permutation_creator, seqs))
    for prod in _product(*perm_creators):
        pmt = sum(prod, ())
        yield pmt


def _product(*iterables, **kwargs):
    if len(iterables) == 0:
        yield ()
    else:
        iterables = iterables * kwargs.get('repeat', 1)
        it = iterables[0]
        for item in it() if callable(it) else iter(it):
            for items in _product(*iterables[1:]):
                yield (item, ) + items


def _permutation_creator(seq):

    def _create():
        return permutations(seq)

    return _create


def _change_indices(bnds, idx_dct):
    if not bnds:
        ret_bnds = frozenset()
    else:
        bkeys, btyps = _bond_keys_and_types(bnds)
        frsts, scnds = _unzip_bond_keys(bkeys)
        ret_frsts = tuple(idx_dct[idx] for idx in frsts)
        ret_scnds = tuple(idx_dct[idx] for idx in scnds)
        ret_bkeys = map(frozenset, zip(ret_frsts, ret_scnds))
        ret_bnds = frozenset(zip(ret_bkeys, btyps))
    return ret_bnds


def _formula(iterable):
    items = tuple(iterable)
    return tuple((item, items.count(item)) for item in sorted(set(items)))


if __name__ == '__main__':
    MGRPH = (('C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
              'H', 'H', 'H', 'H'),
             frozenset([(frozenset([0, 6]), 1), (frozenset([4, 15]), 1),
                        (frozenset([3, 12]), 1), (frozenset([4, 14]), 1),
                        (frozenset([0, 7]), 1), (frozenset([16, 4]), 1),
                        (frozenset([1, 4]), 1), (frozenset([9, 2]), 1),
                        (frozenset([0, 5]), 1), (frozenset([1, 3]), 1),
                        (frozenset([8, 2]), 1), (frozenset([3, 11]), 1),
                        (frozenset([0, 1]), 1), (frozenset([3, 13]), 1),
                        (frozenset([10, 2]), 1), (frozenset([1, 2]), 1)]))
    MGRPH_ = (('C', 'C', 'C', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H',
               'H', 'H', 'H', 'H'),
              frozenset([(frozenset([2, 4]), 1), (frozenset([1, 2]), 1),
                         (frozenset([1, 6]), 1), (frozenset([6, 15]), 1),
                         (frozenset([11, 5]), 1), (frozenset([1, 5]), 1),
                         (frozenset([0, 9]), 1), (frozenset([6, 14]), 1),
                         (frozenset([2, 3]), 1), (frozenset([13, 6]), 1),
                         (frozenset([0, 7]), 1), (frozenset([12, 5]), 1),
                         (frozenset([10, 5]), 1), (frozenset([8, 0]), 1),
                         (frozenset([16, 2]), 1), (frozenset([0, 1]), 1)]))
    print isomorphism(MGRPH, MGRPH_)
