""" parametrized functions for kinetics and thermodynamics
"""
import numpy

GAS_CONSTANT_CAL = 1.9872036


def nasa_enthalpy(temp, cfts):
    """ NASA polynomial enthalpy value in kcal mol-1
    """
    cf1, cf2, cf3, cf4, cf5, cf6, _ = cfts
    enth = GAS_CONSTANT_CAL * 1e3 * (
        cf1 * numpy.power(temp, 1) / 1. +
        cf2 * numpy.power(temp, 2) / 2. +
        cf3 * numpy.power(temp, 3) / 3. +
        cf4 * numpy.power(temp, 4) / 4. +
        cf5 * numpy.power(temp, 5) / 5. + cf6)
    return enth


def arrhenius_log_rate(temp, cfts):
    """ Arrhenius reaction rate log
    """
    arrh_a, arrh_b, arrh_e = cfts
    k = (numpy.log(arrh_a) + numpy.log(temp) * arrh_b
         - arrh_e / temp / GAS_CONSTANT_CAL) / 2.303
    return k
