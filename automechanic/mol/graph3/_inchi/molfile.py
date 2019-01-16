""" MOLFile V3000 string format
"""
import numpy

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
        LINE = (_V3_PFX +
                _ENTRY(key=I_KEY, fmt='d') +      # index
                _SPACE + _ENTRY(key=ORDER_KEY, fmt='d') +  # order
                _SPACE + _ENTRY(key=I1_KEY, fmt='d') +     # atom1 index
                _SPACE + _ENTRY(key=I2_KEY, fmt='d') +     # atom2 index
                _NEWLINE).format


def from_data(atm_keys, bnd_keys, atm_syms, atm_bnd_vlcs, atm_rad_vlcs,
              bnd_ords, is_chi=False, atm_xyzs=None):
    """ MOLFile string from data
    """
    natms = len(atm_keys)
    nbnds = len(bnd_keys)

    key_map = dict(zip(atm_keys, range(1, natms+1)))
    atm_xyzs = numpy.zeros((natms, 3)) if atm_xyzs is None else atm_xyzs

    counts_line = FMT.COUNTS.LINE(
        **{FMT.COUNTS.NA_KEY: natms,
           FMT.COUNTS.NB_KEY: nbnds,
           FMT.COUNTS.CHI_KEY: is_chi})

    atom_block = ''.join((
        FMT.ATOM.LINE(**{FMT.ATOM.I_KEY: key_map[key],
                         FMT.ATOM.S_KEY: sym,
                         FMT.ATOM.X_KEY: x,
                         FMT.ATOM.Y_KEY: y,
                         FMT.ATOM.Z_KEY: z,
                         FMT.ATOM.VAL_KEY: vlc if vlc else -1,
                         FMT.ATOM.MULT_KEY: rad+1})
        for key, sym, (x, y, z), vlc, rad
        in zip(atm_keys, atm_syms, atm_xyzs, atm_bnd_vlcs, atm_rad_vlcs)))

    bond_block = ''.join((
        FMT.BOND.LINE(**{FMT.BOND.I_KEY: i+1,
                         FMT.BOND.ORDER_KEY: ord_,
                         FMT.BOND.I1_KEY: key_map[min(key)],
                         FMT.BOND.I2_KEY: key_map[max(key)]})
        for i, (key, ord_)
        in enumerate(zip(bnd_keys, bnd_ords))))

    mlf = FMT.STRING(**{FMT.COUNTS_KEY: counts_line,
                        FMT.ATOM_KEY: atom_block,
                        FMT.BOND_KEY: bond_block})

    # for recovering the original keys from those used in the molfile
    key_map_inv = dict(map(reversed, key_map.items()))

    return mlf, key_map_inv