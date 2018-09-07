""" NASA coefficient functions
"""

GAS_CONSTANT = 1.9872036e-3


def enthalpy(temp, cfts):
    """ enthalpy value in kcal mol-1
    """
    cf1, cf2, cf3, cf4, cf5, cf6, _ = cfts
    enth = GAS_CONSTANT * (cf1 * temp ** 1 / 1. + cf2 * temp ** 2 / 2. +
                           cf3 * temp ** 3 / 3. + cf4 * temp ** 4 / 4. +
                           cf5 * temp ** 5 / 5. + cf6)
    return enth
