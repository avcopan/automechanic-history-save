""" test the automechanic.table module
"""
import numpy
from automechanic import table2 as table


def test__column_types():
    """ test table.column_types()
    """
    keys = ('a', 'b', 'c')
    typs = tuple((numpy.int32, numpy.float64, numpy.dtype('O')))
    rows = ((1, 1., 'x'), (2, 2., 'y'), (3, 3., 'z'))
    tbl = table.from_rows(rows, keys, typs)
    assert table.column_types(tbl) == typs


def test__column_keys():
    """ test table.column_keys()
    """
    keys = ('a', 'b', 'c')
    typs = tuple((numpy.int32, numpy.float64, numpy.dtype('O')))
    rows = ((1, 1., 'x'), (2, 2., 'y'), (3, 3., 'z'))
    tbl = table.from_rows(rows, keys, typs)
    assert table.column_keys(tbl) == keys


def test__sql_select():
    """ test table.sql_select()
    """
    keys = ('a', 'b', 'c')
    typs = tuple((numpy.int32, numpy.float64, numpy.dtype('O')))
    rows = numpy.array(((1, 1., 'x'), (2, 2., 'y'), (3, 3., 'z')),
                       dtype=numpy.dtype('O'))
    tbl = table.from_rows(rows, keys, typs)
    col_a, col_c = table.table_columns(table.sql_select(tbl, 'a', 'c'))
    assert numpy.array_equiv(col_a, (1, 2, 3))
    assert numpy.array_equiv(col_c, ('x', 'y', 'z'))


if __name__ == '__main__':
    test__column_types()
    test__column_keys()
    test__sql_select()
