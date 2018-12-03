""" common table specifications
"""
from .. import tab

# species table schema
SPC_NAME_KEY = 'name_'
SPC_NAME_TYP = tab.dt_(str)

SPC_MULT_KEY = 'mult_'
SPC_MULT_TYP = tab.dt_(int)

SPC_ID_SMI_KEY = 'smi_'  # smiles
SPC_ID_ICH_KEY = 'ich_'  # inchi
SPC_ID_XYZ_KEY = 'xyz_'  # xyz file path
SPC_ID_TYP = tab.dt_(str)

NASA_C_TYP = tab.dt_(float)
NASA_C_LO_KEYS = ('nasa_lo_1_', 'nasa_lo_2_', 'nasa_lo_3_', 'nasa_lo_4_',
                  'nasa_lo_5_', 'nasa_lo_6_', 'nasa_lo_7_')
NASA_C_HI_KEYS = ('nasa_hi_1_', 'nasa_hi_2_', 'nasa_hi_3_', 'nasa_hi_4_',
                  'nasa_hi_5_', 'nasa_hi_6_', 'nasa_hi_7_')
NASA_T_TYP = tab.dt_(float)
NASA_T_KEYS = ('t_lo_', 't_hi_', 't_c_')

# reaction table schema
RXN_NAME_KEY = 'name_'
RXN_NAME_TYP = tab.dt_(str)

ARRH_TYP = tab.dt_(float)
ARRH_KEYS = ('arrh_a_', 'arrh_b_', 'arrh_e_')
