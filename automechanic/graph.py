""" graph functions
"""
from itertools import chain
from itertools import product
from itertools import permutations
from itertools import combinations
import networkx
from .atom import valence
from .perm import inverse
from .perm import compose


def atoms(mgrph):
    """ the atoms in this molecular graph
    """
    atms, _ = mgrph
    return atms


def indices(mgrph):
    """ indices of the atoms in thie graph
    """
    return tuple(range(len(atoms(mgrph))))


def formula(mgrph):
    """ the molecular formula of this molecular graph
    """
    atms = atoms(mgrph)
    return dict(_formula(atms))


def heavy_atom_indices(mgrph):
    """ heavy atom indices
    """
    idxs = tuple(idx for idx in indices(mgrph)
                 if atom_symbol(mgrph, idx) != 'H')
    return idxs


def hydrogen_atom_indices(mgrph):
    """ heavy atom indices
    """
    idxs = tuple(idx for idx in indices(mgrph)
                 if atom_symbol(mgrph, idx) == 'H')
    return idxs


def backbone_indices(mgrph):
    """ heavy atoms and hydrogen atoms not connected to heavy atoms
    """
    hv_idxs = heavy_atom_indices(mgrph)
    hy_idxs = tuple(idx for idx in hydrogen_atom_indices(mgrph) if
                    all(idx > nidx and atom_symbol(mgrph, nidx) == 'H'
                        for nidx in atom_neighborhood_indices(mgrph, idx)))
    return hv_idxs + hy_idxs


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


def subgraph(mgrph, idxs):
    """ subgraph of some atoms in this molecular graph
    """
    atms = tuple(atom_symbol(mgrph, idx) for idx in idxs)
    bnds = bonds(mgrph)
    sub_bnds = frozenset([(bkey, btyp) for bkey, btyp in bnds
                          if all(bidx in idxs for bidx in bkey)])
    idx_dct = dict((idx, new_idx) for new_idx, idx in enumerate(idxs))
    bnds = _change_indices(sub_bnds, idx_dct)
    return atms, bnds


def heavy_atom_subgraph(mgrph):
    """ the subgraph of the heavy atoms in this molecular graph
    """
    idxs = heavy_atom_indices(mgrph)
    return subgraph(mgrph, idxs)


def skeleton_graph(mgrph):
    """ a graph of the heavy atoms and their hydrogen counts
    """
    idxs = backbone_indices(mgrph)
    atms, bnds = subgraph(mgrph, idxs)
    hcts = (atom_hydrogen_count(mgrph, idx) for idx in idxs)
    vtcs = tuple(zip(atms, hcts))
    return (vtcs, bnds)


def multibond_keys(mgrph):
    """ bond keys for multibonds in the graph
    """
    bdct = bond_orders(mgrph)
    mbnd_bkeys = tuple(key for key, typ in bdct.items() if typ > 1)
    return mbnd_bkeys


def multibond_opening_resonances(mgrph):
    """ X=Y <-> X.-Y. resonances
    """
    mbnd_bkeys = multibond_keys(mgrph)

    res_mgrphs = (mgrph,)
    res_mgrphs += tuple(increment_bond_order(mgrph, idx1, idx2, incr=-1)
                        for idx1, idx2 in mbnd_bkeys)

    return res_mgrphs


def multibond_forming_resonances(mgrph):
    """ X.-Y. <-> X=Y resonances
    """
    mbnd_bkeys = potential_bond_keys(mgrph)

    if len(mbnd_bkeys) == 1:
        bkey, = mbnd_bkeys
        idx1, idx2 = bkey
        res_mgrphs = (mgrph, increment_bond_order(mgrph, idx1, idx2, incr=+1))
    else:
        seeds = tuple(increment_bond_order(mgrph, idx1, idx2, incr=+1)
                      for idx1, idx2 in mbnd_bkeys)
        res_mgrphs = (mgrph,)
        res_mgrphs += tuple(chain(
            *(multibond_forming_resonances(seed) for seed in seeds)))

    return res_mgrphs


def potential_bond_keys(mgrph):
    """ neighboring radical sites of a molecular graph
    """
    ridxs = radical_sites(mgrph)
    return tuple(frozenset([ridx1, ridx2])
                 for ridx1, ridx2 in combinations(ridxs, 2)
                 if ridx2 in atom_neighborhood_indices(mgrph, ridx1))


def radical_sites(mgrph):
    """ radical sites of a molecular graph
    """
    idxs = tuple(idx for idx in indices(mgrph)
                 if atom_free_electrons(mgrph, idx) > 0)
    return idxs


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


def atom_neighborhood(mgrph, idx):
    """ atoms adjacent to this one in a molecule graph
    """
    atms = tuple(atom_symbol(mgrph, nidx)
                 for nidx in atom_neighborhood_indices(mgrph, idx))
    return atms


def atom_hydrogen_indices(mgrph, idx):
    """ indices of the hydrogens attached to this atom
    """
    idxs = tuple(h_idx for h_idx in atom_neighborhood_indices(mgrph, idx)
                 if atom_symbol(mgrph, h_idx) == 'H')
    return idxs


def atom_hydrogen_count(mgrph, idx):
    """ the number of hydrogens attached to this atom
    """
    return len(atom_hydrogen_indices(mgrph, idx))


def atom_bond_count(mgrph, idx):
    """ bond count of an atom in a molecule graph
    """
    atm_bnds = atom_bonds(mgrph, idx)
    _, atm_btyps = _bond_keys_and_types(atm_bnds)
    return sum(atm_btyps)


def atom_symbol(mgrph, idx):
    """ the atomic symbol of an atom by index
    """
    atms = atoms(mgrph)
    return atms[idx]


def atom_free_electrons(mgrph, idx):
    """ number of unbound valence electrons for an atom in a molecular graph
    """
    atms = atoms(mgrph)
    vlnc = valence(atms[idx])
    bcnt = atom_bond_count(mgrph, idx)
    return vlnc - bcnt


def bond_order(mgrph, idx1, idx2):
    """ the order of a bond
    """
    bdct = bond_orders(mgrph)
    return bdct[frozenset([idx1, idx2])]


def permute_atoms(mgrph, pmt):
    """ permute the atoms in a molecular graph
    """
    atms, bnds = mgrph
    ret_atms = tuple(atms[i] for i in pmt)
    inv_pmt = inverse(pmt)
    idx_dct = dict(enumerate(inv_pmt))
    ret_bnds = _change_indices(bnds, idx_dct)
    return ret_atms, ret_bnds


def skeleton_argsort(mgrph):
    """ backbone indices, followed by hydrogen indices for each in turn
    """
    bb_idxs = backbone_indices(mgrph)
    idxs = tuple(chain(
        bb_idxs, *(atom_hydrogen_indices(mgrph, idx) for idx in bb_idxs)))
    return idxs


def skeleton_sort(mgrph):
    """ backbone atoms, followed by hydrogens for each in turn
    """
    pmt = skeleton_argsort(mgrph)
    return permute_atoms(mgrph, pmt)


def union(mgrph1, mgrph2):
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
    mgrph = union(mgrph1, mgrph2)
    mgrph = increment_bond_order(mgrph, idx1, shifted_idx2, incr=order)
    return mgrph


def isomorphic(mgrph1, mgrph2):
    """ are these graphs isomorphic?
    """
    return bool(isomorphism(mgrph1, mgrph2))


def isomorphism(mgrph1, mgrph2):
    """ isomorphism between graphs, if there is one
    """
    iso = None
    if formula(mgrph1) == formula(mgrph2):
        srt1 = skeleton_argsort(mgrph1)
        srt2 = skeleton_argsort(mgrph2)
        srt_mgrph1 = permute_atoms(mgrph1, srt1)
        srt_mgrph2 = permute_atoms(mgrph2, srt2)
        sk_grph1 = skeleton_graph(srt_mgrph1)
        sk_grph2 = skeleton_graph(srt_mgrph2)
        sk_iso = _isomorphism(sk_grph1, sk_grph2)
        if sk_iso:
            srt_iso = tuple(chain(
                sk_iso,
                *(atom_hydrogen_indices(srt_mgrph1, idx) for idx in sk_iso)))
            iso = compose(inverse(srt2), srt_iso, srt1)
    return iso


def forward_abstraction_indices(qh_mgrph, q_mgrph):
    """ find abstraction indices
    """
    idxs = None
    for idx in radical_sites(q_mgrph):
        qh_mgrph_ = bind_atom(q_mgrph, idx, 'H')
        iso = isomorphism(qh_mgrph, qh_mgrph_)
        if iso:
            qh_idx = int(iso[-1])
            q_idx = int(idx)
            idxs = (qh_idx, q_idx)
    return idxs


def addition_indices(x_mgrph, y_mgrph, xy_mgrph):
    """ find addition indices (all kinds)
    """
    idxs = multibond_addition_indices(x_mgrph, y_mgrph, xy_mgrph)

    if not idxs:
        idxs = radical_addition_indices(x_mgrph, y_mgrph, xy_mgrph)

    return idxs


def radical_addition_indices(x_mgrph, y_mgrph, xy_mgrph):
    """ find indices for additions to radical sites
    """
    idxs = None

    x_idxs = radical_sites(x_mgrph)
    y_idxs = radical_sites(y_mgrph)

    for x_idx, y_idx in product(x_idxs, y_idxs):
        xy_mgrph_ = bind(x_mgrph, y_mgrph, x_idx, y_idx)
        for xy_res_mgrph in multibond_forming_resonances(xy_mgrph_):
            iso = isomorphism(xy_mgrph, xy_res_mgrph)
            if iso:
                natms_y = len(atoms(y_mgrph))
                xy_idx_y = int(iso[y_idx])
                xy_idx_x = int(iso[natms_y + x_idx])
                idxs = (x_idx, y_idx, xy_idx_x, xy_idx_y)

    return idxs


def multibond_addition_indices(x_mgrph, y_mgrph, xy_mgrph):
    """ find indices for additions to multiple bonds
    """
    idxs = None

    y_idxs = radical_sites(y_mgrph)
    x_bkeys = multibond_keys(x_mgrph)
    for y_idx, x_bkey in product(y_idxs, x_bkeys):
        for x_idx, x_idx_other in permutations(x_bkey):
            x_open_mgrph = increment_bond_order(
                x_mgrph, x_idx, x_idx_other, incr=-1)
            xy_open_mgrph = bind(y_mgrph, x_open_mgrph, y_idx, x_idx)
            for xy_res_mgrph in multibond_forming_resonances(xy_open_mgrph):
                iso = isomorphism(xy_mgrph, xy_res_mgrph)
                if iso:
                    natms_y = len(atoms(y_mgrph))
                    xy_idx_y = int(iso[y_idx])
                    xy_idx_x = int(iso[natms_y + x_idx])
                    idxs = (x_idx, y_idx, xy_idx_x, xy_idx_y)

    return idxs


def migration_indices(r_mgrph, p_mgrph):
    """ find migration indices
    """
    idxs = None

    r_idxs = radical_sites(r_mgrph)
    p_idxs = radical_sites(p_mgrph)

    for r_idx, p_idx in product(r_idxs, p_idxs):
        r_mgrph_ = bind_atom(r_mgrph, r_idx, 'H')
        p_mgrph_ = bind_atom(p_mgrph, p_idx, 'H')
        iso = isomorphism(r_mgrph_, p_mgrph_)
        if iso:
            r_idx_h = int(iso[-1])
            r_idx_a = int(r_idx)
            p_idx_h = int(inverse(iso)[-1])
            p_idx_a = int(p_idx)
            idxs = (r_idx_h, r_idx_a, p_idx_h, p_idx_a)

    return idxs


def _other_vertex(bkey, idx):
    ret_idx, = bkey - frozenset([idx])
    return ret_idx


def _bond_keys_and_types(bnds):
    bkeys, btyps = zip(*bnds) if bnds else ((), ())
    return bkeys, btyps


def _unzip_bond_keys(bkeys):
    frsts, scnds = zip(*bkeys) if bkeys else ((), ())
    return frsts, scnds


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


def _isomorphism(mgrph1, mgrph2):
    """ full graph isomorphism
    """
    vdct1 = dict(enumerate(atoms(mgrph1)))
    vdct2 = dict(enumerate(atoms(mgrph2)))
    edct1 = bond_orders(mgrph1)
    edct2 = bond_orders(mgrph2)

    nxgrph1 = _nx_graph(vdct1, edct1)
    nxgrph2 = _nx_graph(vdct2, edct2)
    return _nx_graph_isomorphism(nxgrph1, nxgrph2)


def _nx_graph(vdct, edct):
    nxgrph = networkx.Graph()
    for vkey, vtyp in vdct.items():
        nxgrph.add_node(vkey, typ=vtyp)
    for ekey, etyp in edct.items():
        nxgrph.add_edge(*ekey, typ=etyp)
    return nxgrph


def _nx_graph_isomorphism(nxgrph1, nxgrph2):

    def _same_type(dct1, dct2):
        return dct1['typ'] == dct2['typ']

    matcher = networkx.algorithms.isomorphism.GraphMatcher(
        nxgrph1, nxgrph2, node_match=_same_type, edge_match=_same_type)

    iso = None
    if matcher.is_isomorphic():
        iso_dct = matcher.mapping
        iso_inv = tuple(map(iso_dct.__getitem__, range(len(iso_dct))))
        iso = inverse(iso_inv)

    return iso


if __name__ == '__main__':
    MGRPH = (('C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'),
             frozenset([(frozenset([0, 3]), 1), (frozenset([1, 5]), 1),
                        (frozenset([1, 2]), 1), (frozenset([0, 4]), 1),
                        (frozenset([8, 2]), 1), (frozenset([2, 6]), 1),
                        (frozenset([0, 1]), 2), (frozenset([2, 7]), 1)]))
    PMT_MGRPH = (('H', 'H', 'C', 'H', 'C', 'C', 'H', 'H', 'H'),
                 frozenset([(frozenset([5, 6]), 1), (frozenset([4, 7]), 1),
                            (frozenset([8, 4]), 1), (frozenset([0, 5]), 1),
                            (frozenset([1, 5]), 1), (frozenset([2, 5]), 1),
                            (frozenset([2, 3]), 1), (frozenset([2, 4]), 2)]))
    HATM_MGRPH = (('H',), frozenset([]))
    HMOL_MGRPH = (('H', 'H'), frozenset([(frozenset([0, 1]), 1)]))
    MGRPH1 = union(union(HATM_MGRPH, HMOL_MGRPH), MGRPH)
    MGRPH2 = union(union(PMT_MGRPH, HATM_MGRPH), HMOL_MGRPH)
    ISO = _isomorphism(MGRPH1, MGRPH2)
    print permute_atoms(MGRPH1, ISO) == MGRPH2
    ISO = isomorphism(MGRPH1, MGRPH2)
    print permute_atoms(MGRPH1, ISO) == MGRPH2
