""" functions for working with pandas DataFrames
"""
import pandas


def merge(dfrs, col_key):
    """ merge several dataframes into one, by keys
    """
    print dfrs[0][col_key]


if __name__ == '__main__':
    import numpy
    D1 = pandas.DataFrame({'a': numpy.arange(5),
                           'b': 2 * numpy.ones(5),
                           'c': 11 * numpy.ones(5)})
    D2 = pandas.DataFrame({'a': numpy.arange(2, 7),
                           'b': 3 * numpy.ones(5),
                           'd': 12 * numpy.ones(5)})
    print D1
    print D2
    print pandas.DataFrame.merge(D1, D2, on='a')
