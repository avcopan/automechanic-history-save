""" interface to pandas DataFrames modeled on SQL syntax
"""
# import operator
import numpy
from pandas import Series
from pandas import Int64Index
from pandas import RangeIndex
from pandas import concat as _concat

# types
STR_DTYPE = numpy.dtype('O')
INT_DTYPE = numpy.dtype('int64')
FLOAT_DTYPE = numpy.dtype('float64')
BOOL_DTYPE = numpy.dtype('bool')

# values
EMPTY = numpy.nan


def from_data(data, col_keys, col_typs, row_idxs=None):
    """ create a table from a series of rows
    """
    assert numpy.ndim(data) == 2
    data = numpy.array(data)
    nrows, ncols = numpy.shape(data)
    assert len(col_keys) == len(col_typs) == ncols
    row_idxs = _row_indices(nrows, idxs=row_idxs)
    data_by_columns = numpy.transpose(data)
    cols = [_column_series(col_data, key, typ, row_idxs)
            for col_data, key, typ in zip(data_by_columns, col_keys, col_typs)]
    return _concatenate_column_series(cols)


def column_types(tbl):
    """ table column types
    """
    return tuple(tbl.dtypes)


def data_array(tbl):
    """ table data array
    """
    return tbl.data


def column_keys(tbl):
    """ table columns keys
    """
    return tuple(tbl.columns)


def row_indices(tbl):
    """ table row_indices
    """
    return tuple(tbl.index)


def update_column_by_index(tbl, row_idxs, col_key, vals):
    """ update table column by row index
    """
    tbl.loc[list(row_idxs), col_key] = vals
    return tbl


def sql_select_one(tbl, col_key):
    """ select one column from a table
    """
    return tbl[col_key]


def sql_where_eq(tbl, col_key, val):
    """ select rows from a table with column entries of a given value
    """
    assert col_key in column_keys(tbl)
    return tbl[tbl[col_key] == val]


def sql_where_in(tbl, col_key, vals):
    """ select rows from a table with column entries between two values
    """
    assert col_key in column_keys(tbl)
    return tbl[tbl[col_key].isin(vals)]


def _row_indices(nrows, idxs=None):
    row_idxs = Int64Index(idxs) if idxs is not None else RangeIndex(stop=nrows)
    assert len(row_idxs) == nrows
    return row_idxs


def _column_series(col, key, typ, row_idxs):
    return Series(data=col, name=key, dtype=typ, index=row_idxs)


def _concatenate_column_series(col_sers):
    return _concat(col_sers, axis=1)
