""" functions for working with pandas DataFrames
"""
from itertools import chain
from collections import OrderedDict
from more_itertools import unique_everseen
import pandas


def is_empty_value(val):
    """ check if a table value is missing, by value
    """
    return pandas.isna(val)


def column(table_df, col_key):
    """ table column
    """
    assert col_key in column_keys(table_df)
    return tuple(table_df[col_key])


def columns(table_df, col_keys=None):
    """ columns of a table, in order
    """
    if col_keys is None:
        col_keys = column_keys(table_df)

    cols = tuple(column(table_df, col_key) for col_key in col_keys)
    return cols


def rows(table_df):
    """ rows of a table
    """
    return tuple(iterate_rows(table_df))


def iterate_rows(table_df):
    """ iterate over rows of a table, as tuples to preserve order
    """
    for _, row in table_df.iterrows():
        yield dict(row.items())


def update_column_keys(table_df, col_keys):
    """ update column keys, adding new ones to the end
    """
    col_keys_in = column_keys(table_df)
    col_keys_upd = tuple(col_key for col_key in col_keys
                         if col_key not in col_keys_in)
    col_keys_out = col_keys_in + col_keys_upd
    return table_df.reindex(columns=col_keys_out)


def columns_like(table_df):
    """ return an empty table with columns like this one
    """
    return empty(column_keys(table_df))


def column_keys(table_df):
    """ get the column keys of a table
    """
    return tuple(table_df.columns)


def empty(col_keys):
    """ construct an empty pandas.DataFrame
    """
    return pandas.DataFrame(columns=col_keys)


def from_columns(cols, col_keys):
    """ construct a pandas.DataFrame from columns
    """
    assert len(cols) == len(col_keys)
    col_dct = dict(zip(col_keys, cols))
    return pandas.DataFrame(col_dct, columns=col_keys)


def from_rows(row_dcts, col_keys):
    """ construct a pandas.DataFrame from rows
    """
    return pandas.DataFrame(list(row_dcts), columns=col_keys)


def append_rows(table_df, row_dcts):
    """ append rows to a table
    """
    col_keys = column_keys(table_df)
    row_dcts = tuple(rows(table_df)) + tuple(row_dcts)
    return from_rows(row_dcts, col_keys=col_keys)


def append_columns(table_df, cols, col_keys):
    """ append columns to a table
    """
    col_keys = tuple(column_keys(table_df)) + tuple(col_keys)
    cols = tuple(columns(table_df)) + tuple(cols)
    return from_columns(cols, col_keys=col_keys)


def reindex(table_df):
    """ add/overwrite 'index' column with a range index
    """
    table_df['index'] = range(len(table_df))
    table_df = move_column_to_front(table_df, 'index')
    return table_df


def sort(table_df, col_key, descending=False):
    """ sort a table by column
    """
    table_df = pandas.DataFrame.sort_values(
        table_df, by=col_key, ascending=not descending)
    return table_df


def intersect(table_dfs, col_key):
    """ intsc tables by column
    """
    col_key_vals = list(unique_everseen(chain(*(
        table_df[col_key] for table_df in table_dfs))))
    lookup_dcts = [lookup_dictionary(table_df, col_key)
                   for table_df in table_dfs]

    intscd_rows = []
    for val in col_key_vals:
        row = {}
        if val and all(val in lookup_dct for lookup_dct in lookup_dcts):
            for lookup_dct in lookup_dcts:
                row.update(lookup_dct[val])
            intscd_rows.append(row)
    intscd_col_keys = list(unique_everseen(chain(*table_dfs)))
    intscd_df = pandas.DataFrame.from_dict(intscd_rows)[intscd_col_keys]
    return intscd_df


def merge(table_dfs, col_key):
    """ merge tables by column
    """
    col_key_vals = list(unique_everseen(chain(*(
        table_df[col_key] for table_df in table_dfs))))
    lookup_dcts = [lookup_dictionary(table_df, col_key)
                   for table_df in table_dfs]

    merged_rows = []
    for val in col_key_vals:
        row = {col_key: val}
        for lookup_dct in lookup_dcts:
            if val in lookup_dct:
                row.update(lookup_dct[val])
        merged_rows.append(row)
    merged_col_keys = list(unique_everseen(chain(*table_dfs)))
    merged_df = pandas.DataFrame.from_dict(merged_rows)[merged_col_keys]
    return merged_df


def lookup_dictionary(table_df, col_key):
    """ look up table rows by column value
    """
    lkp_col_vals = column(table_df, col_key)
    rows_ = rows(table_df)
    return OrderedDict(zip(lkp_col_vals, rows_))


def from_lookup_dictionary(table_lkp, col_keys):
    """ rectonstruct a table from a lookup dictionary
    """
    rows_ = table_lkp.values()
    return from_rows(rows_, col_keys=col_keys)


def lookup_update(table_df, lookup_item, update_item):
    """ lookup row by column value and update an entry
    """
    col_keys = column_keys(table_df)
    lkp_col_key, lkp_col_val = lookup_item
    upd_col_key, upd_col_val = update_item
    assert lkp_col_key in col_keys and upd_col_key in col_keys
    table_lkp = lookup_dictionary(table_df, lkp_col_key)
    table_lkp[lkp_col_val][upd_col_key] = upd_col_val
    return from_lookup_dictionary(table_lkp, col_keys)


def move_column_to_front(table_df, col_key):
    """ move a column to a new position
    """
    assert col_key in table_df
    col_keys = [col_key] + list(key for key in table_df.columns
                                if key != col_key)
    return table_df[col_keys]


if __name__ == '__main__':
    TAB = empty(('a', 'b'))
    TAB = append_rows(TAB, [{'a': 1, 'b': 2}])
    TAB = append_rows(TAB, [{'a': 3, 'b': 4}])
    print TAB
