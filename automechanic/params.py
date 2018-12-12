""" automechanic run parameters
"""
from . import tab

SPC_FS_DIR_NAME = 'SPC'

# species table schema
SPC_NAME_KEY = 'name'
SPC_NAME_TYP = tab.dt_(str)

SPC_MULT_KEY = 'mult'
SPC_MULT_TYP = tab.dt_(int)

SPC_ID_SMI_KEY = 'smiles'  # smiles
SPC_ID_ICH_KEY = 'inchi'  # inchi
SPC_ID_XYZ_KEY = 'xyz'  # xyz file path
SPC_ID_TYP = tab.dt_(str)

SPC_FS_PATH_KEY = 'fs_path'
SPC_FS_PATH_TYP = tab.dt_(str)

NASA_C_TYP = tab.dt_(float)
NASA_C_LO_KEYS = ('nasa_lo_1', 'nasa_lo_2', 'nasa_lo_3', 'nasa_lo_4',
                  'nasa_lo_5', 'nasa_lo_6', 'nasa_lo_7')
NASA_C_HI_KEYS = ('nasa_hi_1', 'nasa_hi_2', 'nasa_hi_3', 'nasa_hi_4',
                  'nasa_hi_5', 'nasa_hi_6', 'nasa_hi_7')
NASA_T_TYP = tab.dt_(float)
NASA_T_KEYS = ('t_lo', 't_hi', 't_c')

# species keys and options
PICK_STEREO = 'pick'
EXPAND_STEREO = 'expand'

# reaction table schema
RXN_NAME_KEY = 'name'
RXN_NAME_TYP = tab.dt_(str)

RXN_FS_PATH_KEY = 'fs_path'
RXN_FS_PATH_TYP = tab.dt_(str)

ARRH_TYP = tab.dt_(float)
ARRH_KEYS = ('arrh_a', 'arrh_b', 'arrh_e')
