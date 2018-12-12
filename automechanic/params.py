""" automechanic run parameters
"""
from . import tab

_NAME_KEY = 'name'
_NAME_TYP = tab.dt_(str)
_ID_ICH_KEY = 'inchi'
_ID_SMI_KEY = 'smiles'
_ID_TYP = tab.dt_(str)
_MULT_KEY = 'mult'
_MULT_TYP = tab.dt_(int)
_FILESYSTEM_PATH_KEY = 'path'
_FILESYSTEM_PATH_TYP = tab.dt_(str)


class SPC():
    """ species parameters
    """
    ID_ICH_KEY = _ID_ICH_KEY
    ID_SMI_KEY = _ID_SMI_KEY
    MULT_KEY = _MULT_KEY

    FILESYSTEM_DIR_NAME = 'SPC'

    PICK_STEREO = 'pick'
    EXPAND_STEREO = 'expand'

    class TAB():
        """ species table parameters
        """
        NAME_KEY = _NAME_KEY
        NAME_TYP = _NAME_TYP

        ID_TYP = _ID_TYP

        MULT_TYP = _MULT_TYP

        FILESYSTEM_PATH_KEY = _FILESYSTEM_PATH_KEY
        FILESYSTEM_PATH_TYP = _FILESYSTEM_PATH_TYP

        NASA_C_TYP = tab.dt_(float)
        NASA_C_LO_KEYS = ('nasa_lo_1', 'nasa_lo_2', 'nasa_lo_3', 'nasa_lo_4',
                          'nasa_lo_5', 'nasa_lo_6', 'nasa_lo_7')
        NASA_C_HI_KEYS = ('nasa_hi_1', 'nasa_hi_2', 'nasa_hi_3', 'nasa_hi_4',
                          'nasa_hi_5', 'nasa_hi_6', 'nasa_hi_7')
        NASA_T_TYP = tab.dt_(float)
        NASA_T_KEYS = ('t_lo', 't_hi', 't_c')


class RXN():
    """ reaction parameters
    """
    MULT_KEY = _MULT_KEY
    ID_ICH_KEY = _ID_ICH_KEY

    FILESYSTEM_DIR_NAME = 'RXN'

    class TAB():
        """ species table parameters
        """
        NAME_KEY = _NAME_KEY
        NAME_TYP = _NAME_TYP

        ID_TYP = _ID_TYP

        MULT_TYP = _MULT_TYP

        FILESYSTEM_PATH_KEY = _FILESYSTEM_PATH_KEY
        FILESYSTEM_PATH_TYP = _FILESYSTEM_PATH_TYP

        ARRH_TYP = tab.dt_(float)
        ARRH_KEYS = ('arrh_a', 'arrh_b', 'arrh_e')
