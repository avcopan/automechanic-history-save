""" interface to pandas for tabular data
"""
from collections import OrderedDict as _OrderedDict
from itertools import starmap as _starmap
import numpy
from pandas import Series as _Series
from pandas import read_csv as _read_csv
from pandas import DataFrame as _DataFrame
from pandas import Int64Index as _Int64Index
from pandas import RangeIndex as _RangeIndex

IDX_KEY = '_i'


# table creation routines
def from_records(vals, keys, typs=None, idxs=None):
    """ create a table from a series of records
    """
    assert numpy.ndim(vals) == 2
    nrows, ncols = numpy.shape(vals)
    idxs = (_RangeIndex(stop=nrows, name=IDX_KEY) if idxs is None
            else _Int64Index(idxs, name=IDX_KEY))
    typs = tuple(map(dt_, (object,) * ncols if typs is None else typs))
    assert len(keys) == len(typs) == ncols and len(idxs) == nrows
    assert IDX_KEY not in keys

    cols = numpy.transpose(vals)
    data = _OrderedDict((key, _Series(data=col, dtype=typ, index=idxs))
                        for col, key, typ in zip(cols, keys, typs))

    return _DataFrame(data=data, index=idxs)


# def from_fancy_records(fancy_vals, fancy_keys, typs=None, idxs=None):
#     """
#     """
#     pass


def from_starmap(tbl, func, iter_keys, keys, typs=None):
    """ starmap columns and create a table from the resulting values
    """
    vals = list(
        _starmap(func, _iterate_fancy_select(tbl, *iter_keys)))
    return from_records(vals, keys, typs=typs, idxs=idxs_(tbl))


# table I/O
def read_csv(file_path):
    """ read table from a CSV file
    """
    tbl = _read_csv(file_path)
    if IDX_KEY in keys_(tbl):
        tbl.set_index(IDX_KEY, inplace=True)
    else:
        tbl.index.rename(IDX_KEY, inplace=True)
    return tbl


def write_csv(tbl, file_path):
    """ write table to a CSV file
    """
    tbl.to_csv(file_path)


# table modification routines
def update(tbl1, tbl2):
    """ write the values in table 2 over those in table 1

    indices of table 2 must be a subset of those in table 1, but it may
    contain new columns to be added to table 1
    """
    assert set(idxs_(tbl1)) >= set(idxs_(tbl2))
    tbl = tbl1.copy(deep=True)
    for key in tbl2:
        tbl.loc[list(idxs_(tbl2)), key] = tbl2[key]
    return tbl


# data types
def dt_(typ):
    """ return data type conversion
    """
    np_typ = numpy.dtype(typ)
    if numpy.issubdtype(np_typ, numpy.bool_):
        typ = numpy.dtype('bool')
    elif numpy.issubdtype(np_typ, numpy.integer):
        typ = numpy.dtype('int64')
    elif numpy.issubdtype(np_typ, numpy.number):
        typ = numpy.dtype('float64')
    else:
        typ = numpy.dtype('O')
    return typ


# table properties
def idxs_(tbl):
    """ table indices
    """
    return tuple(tbl.index)


def keys_(tbl):
    """ table keys
    """
    return tuple(tbl.columns)


def vals_(tbl):
    """ table values
    """
    return numpy.array(tbl.values)


def typs_(tbl):
    """ table column types
    """
    return tuple(tbl.dtypes)


# SQL functions
def select(tbl, *keys):
    """ one or more columns from a table, ordered
    """
    return tbl.__getitem__(list(keys))


def where(tbl, bool_):
    """ one or more rows from a table, unordered
    """
    return tbl.__getitem__(bool_)


# boolean operators to go along with the `where` function
def lt_(tbl, key, val):
    """ less than
    """
    return tbl[key].lt(val)


def eq_(tbl, key, val):
    """ equals
    """
    return tbl[key].eq(val)


def gt_(tbl, key, val):
    """ greater than
    """
    return tbl[key].gt(val)


def isin_(tbl, key, val):
    """ is in
    """
    return tbl[key].isin(val)


def not_(bool_):
    """ not
    """
    return ~bool_


def and_(bool1, bool2):
    """ and
    """
    return bool1 & bool2


def or_(bool1, bool2):
    """ or
    """
    return bool1 | bool2


# iterators
def iter_(tbl):
    """ yield row values
    """
    for vals in tbl.itertuples(index=False, name=None):
        yield vals


def enum_(tbl):
    """ yield row values with their indices
    """
    for tup in tbl.itertuples(index=True, name=None):
        idx, vals = tup[0], tup[1:]
        yield (idx, vals)


def _iterate_fancy_select(tbl, *keys_and_groups):
    """ yield fancily-selected row values
    """
    iterators = map(_iter_values, _fancy_select(tbl, *keys_and_groups))
    for fancy_vals in zip(*iterators):
        yield fancy_vals


def _fancy_select(tbl, *keys_and_groups):
    """ select individual colums and groups of columns
    """
    return tuple(select(tbl, slct)[slct] if isinstance(slct, str)
                 else select(tbl, *slct) for slct in keys_and_groups)


def _iter_values(tbl):
    """ iterate over values in a table or series
    """
    if hasattr(tbl, 'itertuples'):
        for vals in tbl.itertuples(index=False, name=None):
            yield vals
    elif hasattr(tbl, 'iteritems'):
        for _, val in tbl.iteritems():
            yield val
