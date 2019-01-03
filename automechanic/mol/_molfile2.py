""" MOLFile V3000 string format
"""
import numpy
from .graph2._shared import atom_keys
from .graph2._shared import bond_keys
from .graph2._shared import atom_symbols
from .graph2._shared import atom_neighbor_keys
from .graph2._shared import relabel
from .graph2._shared import highspin_resonance_graph
from .graph2._perm import parity as _parity
from .graph2.res import atom_bond_valences as _rgr_atom_bond_valences
from .graph2.res import atom_radical_valences as _rgr_atom_radical_valences
from .graph2.res import bond_orders as _rgr_bond_orders
from .graph2.stereo import atom_parities as _sgr_atom_parities
from .graph2.stereo import is_chiral as _sgr_is_chiral
from .graph2._dict import values_by_key as _values_by_key
from .atom import nuclear_charge as _atom_nuclear_charge

_V3_PFX = 'M  V30 '
_SPACE = ' '
_ENTRY = '{{{key}:{fmt}}}'.format
_ZERO = '0'
_NEWLINE = '\n'


class FMT():
    """ MOLFile V3000 format specifications """
    _HEAD = ('' + _NEWLINE +
             '  automech  2D grid' + _NEWLINE +
             '' + _NEWLINE +
             '  0  0  0  0  0  0  0  0  0  0999 V3000' + _NEWLINE)
    _FOOT = 'M  END' + _NEWLINE
    _CTAB = 'CTAB'
    _ATOM = 'ATOM'
    _BOND = 'BOND'
    _BEGIN = (_V3_PFX + 'BEGIN {:s}' + _NEWLINE).format
    _END = (_V3_PFX + 'END {:s}' + _NEWLINE).format
    COUNTS_KEY = 'counts'
    ATOM_KEY = 'atom'
    BOND_KEY = 'bond'
    STRING = (_HEAD + _BEGIN(_CTAB) +
              _ENTRY(key=COUNTS_KEY, fmt='s') +
              _BEGIN(_ATOM) + _ENTRY(key=ATOM_KEY, fmt='s') + _END(_ATOM) +
              _BEGIN(_BOND) + _ENTRY(key=BOND_KEY, fmt='s') + _END(_BOND) +
              _END(_CTAB) + _FOOT).format

    class COUNTS():
        """ _ """
        NA_KEY = 'na'
        NB_KEY = 'nb'
        CHI_KEY = 'chiral'
        LINE = (_V3_PFX + 'COUNTS' +
                _SPACE + _ENTRY(key=NA_KEY, fmt='d') +   # atom count
                _SPACE + _ENTRY(key=NB_KEY, fmt='d') +   # bond count
                _SPACE + _ZERO +                         # no Sgroups
                _SPACE + _ZERO +                         # no 3d constraints
                _SPACE + _ENTRY(key=CHI_KEY, fmt='d') +  # is chiral?
                _NEWLINE).format

    class ATOM():
        """ _ """
        I_KEY = 'i'
        S_KEY = 's'
        X_KEY = 'x'
        Y_KEY = 'y'
        Z_KEY = 'z'
        MULT_KEY = 'mult'
        VAL_KEY = 'valence'
        LINE = (_V3_PFX +
                _ENTRY(key=I_KEY, fmt='d') +   # index
                _SPACE + _ENTRY(key=S_KEY, fmt='s') +   # symbol
                _SPACE + _ENTRY(key=X_KEY, fmt='.3f') +     # x coordinate
                _SPACE + _ENTRY(key=Y_KEY, fmt='.3f') +     # y coordinate
                _SPACE + _ENTRY(key=Z_KEY, fmt='.3f') +     # z coordinate
                _SPACE + 'RAD=' + _ENTRY(key=MULT_KEY, fmt='d') +
                _SPACE + 'VAL=' + _ENTRY(key=VAL_KEY, fmt='d') +
                _NEWLINE).format

    class BOND():
        """ _ """
        I_KEY = 'i'
        ORDER_KEY = 'order'
        I1_KEY = 'i1'
        I2_KEY = 'i2'
        CFG_KEY = 'stereo_config'
        LINE = (_V3_PFX +
                _ENTRY(key=I_KEY, fmt='d') +      # index
                _SPACE + _ENTRY(key=ORDER_KEY, fmt='d') +  # order
                _SPACE + _ENTRY(key=I1_KEY, fmt='d') +     # atom1 index
                _SPACE + _ENTRY(key=I2_KEY, fmt='d') +     # atom2 index
                _SPACE + 'CFG=' + _ENTRY(key=CFG_KEY, fmt='d') +
                _NEWLINE).format


def from_stereo_graph(sgr):
    """ InChI string from a stereo graph
    """
    sgr = _relabel_with_one_indexing(sgr)
    atm_keys = atom_keys(sgr)
    bnd_keys = bond_keys(sgr)

    # counts line
    natms = len(atm_keys)
    nbnds = len(bnd_keys)
    is_chi = _sgr_is_chiral(sgr)
    counts_line = FMT.COUNTS.LINE(
        **{FMT.COUNTS.NA_KEY: natms,
           FMT.COUNTS.NB_KEY: nbnds,
           FMT.COUNTS.CHI_KEY: is_chi})

    # atom block
    rgr = highspin_resonance_graph(sgr)
    atm_bnd_vlcs = _values_by_key(_rgr_atom_bond_valences(rgr), atm_keys)
    atm_rad_vlcs = _values_by_key(_rgr_atom_radical_valences(rgr), atm_keys,
                                  fill=0)

    atm_syms = _values_by_key(atom_symbols(sgr), atm_keys)
    atm_mlf_vlcs = [atm_bnd_vlc if atm_bnd_vlc != 0 else -1
                    for atm_bnd_vlc in atm_bnd_vlcs]
    atm_mlf_mults = numpy.add(atm_rad_vlcs, 1)
    atm_xyzs = numpy.zeros((natms, 3))

    atom_block = ''.join((
        FMT.ATOM.LINE(**{FMT.ATOM.I_KEY: key+1,
                         FMT.ATOM.S_KEY: sym,
                         FMT.ATOM.X_KEY: x,
                         FMT.ATOM.Y_KEY: y,
                         FMT.ATOM.Z_KEY: z,
                         FMT.ATOM.VAL_KEY: vlc,
                         FMT.ATOM.MULT_KEY: mult})
        for key, sym, (x, y, z), vlc, mult
        in zip(atm_keys, atm_syms, atm_xyzs, atm_mlf_vlcs, atm_mlf_mults)))

    # bond block
    bnd_mlf_cfgs = _values_by_key(_molfile_bond_configurations(sgr), bnd_keys,
                                  fill=0)
    bnd_ords = _values_by_key(_rgr_bond_orders(rgr), bnd_keys)
    bond_block = ''.join((
        FMT.BOND.LINE(**{FMT.BOND.I_KEY: i+1,
                         FMT.BOND.ORDER_KEY: ord_,
                         FMT.BOND.I1_KEY: min(key)+1,
                         FMT.BOND.I2_KEY: max(key)+1,
                         FMT.BOND.CFG_KEY: cfg})
        for i, (key, ord_, cfg)
        in enumerate(zip(bnd_keys, bnd_ords, bnd_mlf_cfgs))))

    mlf = FMT.STRING(**{FMT.COUNTS_KEY: counts_line,
                        FMT.ATOM_KEY: atom_block,
                        FMT.BOND_KEY: bond_block})
    return mlf


def _relabel_with_one_indexing(sgr):
    atm_keys = atom_keys(sgr)
    new_atm_keys = range(1, len(atm_keys)+1)
    return relabel(sgr, dict(zip(atm_keys, new_atm_keys)))


def _molfile_bond_configurations(sgr):
    atm_par_dct = _sgr_atom_parities(sgr)
    atm_keys = atm_par_dct.keys()
    # atm_pars = atm_par_dct.values()
    # atm_cfgs = [3 if atm_par is True else 1 for atm_par in atm_pars]
    atm_ngb_keys_lst = _values_by_key(atom_neighbor_keys(sgr), atm_keys)
    atm_ngb_key_lst = [
        atm_ngb_keys[-2] if len(atm_ngb_keys) == 4 else
        atm_ngb_keys[-1] if len(atm_ngb_keys) == 3 else
        None
        for atm_ngb_keys in atm_ngb_keys_lst]

    print(atm_ngb_keys_lst)
    print(atm_ngb_key_lst)

    atm_sym_dct = atom_symbols(sgr)
    atm_ngb_chgs_lst = tuple(
        tuple(map(_atom_nuclear_charge,
                  _values_by_key(atm_sym_dct, atm_ngb_keys)))
        for atm_ngb_keys in atm_ngb_keys_lst)
    atm_ngb_srts = tuple(
        numpy.argsort(atm_ngb_chgs)
        for atm_ngb_keys, atm_ngb_chgs
        in zip(atm_ngb_keys_lst, atm_ngb_chgs_lst))
    atm_cfgs = [3 if _parity(atm_ngb_srt, sorted(atm_ngb_srt)) else 1
                for atm_ngb_srt in atm_ngb_srts]
    print(atm_ngb_srts)
    print(atm_cfgs)

    print(atm_ngb_key_lst)
    return {frozenset({atm_key, atm_ngb_key}): atm_cfg
            for atm_key, atm_ngb_key, atm_cfg
            in zip(atm_keys, atm_ngb_key_lst, atm_cfgs)}


# def _molfile_bond_configurations(sgr):
#     atm_keys, atm_pars = zip(*_sgr_atom_parities(sgr).items())
#     atm_cfgs = [1 if atm_par is True else 3 for atm_par in atm_pars]
#     print(atm_cfgs)
#     atm_ngb_keys_dct = atom_neighbor_keys(sgr)
#     atm_ngb_keys = tuple(map(next, map(iter, map(
#         sorted, _values_by_key(atm_ngb_keys_dct, atm_keys)))))
#     return {frozenset({atm_key, atm_ngb_key}): atm_cfg
#             for atm_key, atm_ngb_key, atm_cfg
#             in zip(atm_keys, atm_ngb_keys, atm_cfgs)}
