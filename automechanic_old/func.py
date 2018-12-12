""" parametrized functions for kinetics and thermodynamics
"""
import numpy

GAS_CONSTANT_CAL = 1.9872036


def nasa_enthalpy_function(cfts_lo, cfts_hi, temp_spec):
    """ enthalpy in kcal mol-1
    """
    enth_over_rt_ = nasa_enthalpy_polynomial(cfts_lo, cfts_hi, temp_spec)

    def _enthalpy(temp_):
        rt_ = GAS_CONSTANT_CAL * 1e3 * temp_
        return enth_over_rt_(temp_) * rt_

    return _enthalpy


def nasa_gibbs_polynomial(cfts_lo, cfts_hi, temp_spec):
    """ callable NASA Gibbs polynomial
    """
    enth_over_rt_ = nasa_enthalpy_polynomial(cfts_lo, cfts_hi, temp_spec)
    entr_over_r_ = nasa_entropy_polynomial(cfts_lo, cfts_hi, temp_spec)

    def _gibbs_over_rt(temp_):
        return enth_over_rt_(temp_) - entr_over_r_(temp_)

    return _gibbs_over_rt


def nasa_enthalpy_polynomial(cfts_lo, cfts_hi, temp_spec):
    """ callable NASA enthalpy polynomial
    """
    return _vectorized_nasa_function(
        nasa_enthalpy_polynomial_value, cfts_lo, cfts_hi, temp_spec)


def nasa_entropy_polynomial(cfts_lo, cfts_hi, temp_spec):
    """ callable NASA entropy polynomial
    """
    return _vectorized_nasa_function(
        nasa_entropy_polynomial_value, cfts_lo, cfts_hi, temp_spec)


def nasa_enthalpy_polynomial_value(temp, cfts):
    """ NASA polynomial enthalpy value
    """
    cf1, cf2, cf3, cf4, cf5, cf6, _ = cfts
    enth_over_rt = (cf1 * numpy.power(temp, 0) / 1. +
                    cf2 * numpy.power(temp, 1) / 2. +
                    cf3 * numpy.power(temp, 2) / 3. +
                    cf4 * numpy.power(temp, 3) / 4. +
                    cf5 * numpy.power(temp, 4) / 5. +
                    numpy.divide(cf6, temp))
    return enth_over_rt


def nasa_entropy_polynomial_value(temp, cfts):
    """ NASA polynomial entropy value
    """
    cf1, cf2, cf3, cf4, cf5, _, cf7 = cfts
    entr_over_r = (cf1 * numpy.log(temp) +
                   cf2 * numpy.power(temp, 1) / 1. +
                   cf3 * numpy.power(temp, 2) / 2. +
                   cf4 * numpy.power(temp, 3) / 3. +
                   cf5 * numpy.power(temp, 4) / 4. +
                   cf7)
    return entr_over_r


def arrhenius_log_rate(temp, cfts):
    """ Arrhenius reaction rate log
    """
    arrh_a, arrh_b, arrh_e = cfts
    k = (numpy.log(arrh_a) + numpy.log(temp) * arrh_b
         - arrh_e / temp / GAS_CONSTANT_CAL) / numpy.log(10.)
    return k


# helpers
def _select_nasa_coefficients(temp, cfts_lo, cfts_hi, temp_spec):
    temp_com, _, _ = temp_spec
    cfts = cfts_lo if temp <= temp_com else cfts_hi
    return cfts


def _vectorized_nasa_function(nasa_func, cfts_lo, cfts_hi, temp_spec):

    def _func(temp_):
        cfts = _select_nasa_coefficients(temp_, cfts_lo, cfts_hi, temp_spec)
        return nasa_func(temp_, cfts)

    return numpy.vectorize(_func)
