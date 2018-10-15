""" interface to pandas DataFrames modeled on SQL syntax
"""
import numpy
from pandas import Series
from pandas import concat as _concat


def from_rows(rows, keys, typs):
    """ create a table from a series of rows
    """
    cols = numpy.transpose(numpy.array(rows))
    return from_columns(cols, keys, typs)


def from_columns(cols, keys, typs):
    """ create a table from a series columns
    """
    assert len(cols) == len(keys) == len(typs)
    col_series = [Series(data=col, name=key, dtype=typ)
                  for col, key, typ in zip(cols, keys, typs)]
    return _concatenate_column_series(col_series)


def column_types(tbl):
    """ table column types
    """
    return tuple(tbl.dtypes)


def column_keys(tbl):
    """ table columns keys
    """
    return tuple(tbl.columns)


def table_rows(tbl):
    """ rows of a table
    """
    return tbl.values


def table_columns(tbl):
    """ columns of a table
    """
    return numpy.transpose(table_rows(tbl))


def sql_select(tbl, *keys):
    """ select columns from a table
    """
    return tbl[list(keys)]


def _concatenate_column_series(srys):
    return _concat(srys, axis=1)
