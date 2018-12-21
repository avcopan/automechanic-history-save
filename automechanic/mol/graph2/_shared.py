""" shared molecular graph functions
"""

SYMBOL_POS = 0
HCOUNT_POS = 1
PARITY_POS = 2


def _from_data(atm_keys, atm_vals, bnd_keys, bnd_vals):
    """ construct molecular graph from data
    """
    assert len(atm_keys) == len(atm_vals)
    assert len(bnd_keys) == len(bnd_vals)
    assert all(isinstance(bnd_key, frozenset) and len(bnd_key) == 2
               for bnd_key in bnd_keys)
    assert frozenset.union(*bnd_keys) <= set(atm_keys)
    atms = dict(zip(atm_keys, atm_vals))
    bnds = dict(zip(bnd_keys, bnd_vals))
    xgr = (atms, bnds)
    return xgr


def atoms(xgr):
    """ atoms
    """
    atms, _ = xgr
    return atms


def bonds(xgr):
    """ bonds
    """
    _, bnds = xgr
    return bnds


def atom_keys(xgr):
    """ sorted atom keys
    """
    atm_keys, _ = zip(*_sorted_atoms_list(xgr))
    return atm_keys


def _atom_values(xgr):
    """ atom values, sorted by atom key
    """
    _, atm_vals = zip(*_sorted_atoms_list(xgr))
    return atm_vals


def atom_symbols(xgr):
    """ atomic symbols, sorted by atom key
    """
    atm_syms = list(zip(*_atom_values(xgr)))[SYMBOL_POS]
    return atm_syms


def atom_hydrogen_counts(xgr):
    """ atom hydrogen counts, sorted by atom key
    """
    atm_hyd_cnts = list(zip(*_atom_values(xgr)))[HCOUNT_POS]
    return atm_hyd_cnts


def atom_stereo_parities(xgr):
    """ atom hydrogen counts, sorted by atom key
    """
    atm_pars = list(zip(*_atom_values(xgr)))[PARITY_POS]
    return atm_pars


def bond_keys(xgr):
    """ sorted bond keys
    """
    atm_keys, _ = zip(*_sorted_bonds_list(xgr))
    return atm_keys


def _bond_values(xgr):
    """ bond values, sorted by key
    """
    _, atm_vals = zip(*_sorted_bonds_list(xgr))
    return atm_vals


def _sorted_atoms_list(xgr):
    """ sorted list of atoms
    """
    return sorted(atoms(xgr).items(), key=lambda x: x[0])


def _sorted_bonds_list(xgr):
    """ sorted list of bonds
    """
    return sorted(bonds(xgr).items(), key=lambda x: tuple(sorted(x[0])))
