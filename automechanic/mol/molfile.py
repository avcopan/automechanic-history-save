""" functions operating on MOLFile V3000 strings
"""
from ._irdkit import from_molfile as _rdm_from_molfile
from ._irdkit import to_inchi as _rdm_to_inchi

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
        CFG_KEY = 'stereo_config'
        LINE = (_V3_PFX +
                _ENTRY(key=I_KEY, fmt='d') +   # index
                _SPACE + _ENTRY(key=S_KEY, fmt='s') +   # symbol
                _SPACE + _ENTRY(key=X_KEY, fmt='.3f') +     # x coordinate
                _SPACE + _ENTRY(key=Y_KEY, fmt='.3f') +     # y coordinate
                _SPACE + _ENTRY(key=Z_KEY, fmt='.3f') +     # z coordinate
                _SPACE + 'RAD=' + _ENTRY(key=MULT_KEY, fmt='d') +
                _SPACE + 'VAL=' + _ENTRY(key=VAL_KEY, fmt='d') +
                _SPACE + 'CFG=' + _ENTRY(key=CFG_KEY, fmt='d') +
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


def inchi(mlf):
    """ InChI string from a MOLFile string
    """
    rdm = _rdm_from_molfile(mlf)
    ich = _rdm_to_inchi(rdm)
    return ich
