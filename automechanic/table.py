""" functions for working with pandas DataFrames
"""
from itertools import chain
from more_itertools import unique_everseen
import pandas


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
    other_col_keys = [key for key in table_df.columns
                      if key != col_key]
    _, rows = zip(*table_df[other_col_keys].iterrows())
    lookup_keys = table_df[col_key]
    lookup_values = map(dict, rows)
    return dict(zip(lookup_keys, lookup_values))


def move_column_to_front(table_df, col_key):
    """ move a column to a new position
    """
    assert col_key in table_df
    col_keys = [col_key] + list(key for key in table_df.columns
                                if key != col_key)
    return table_df[col_keys]
