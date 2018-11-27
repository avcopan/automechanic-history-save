""" I/O helpers
"""
from __future__ import unicode_literals
from builtins import open
import os
import time


def timestamp_if_exists(fpath):
    """ open a file, avoiding overwrites if requested
    """
    if os.path.isfile(fpath):
        time_stamp = time.strftime("%Y%m%d-%H%M%S")
        new_fpath = "{:s}_{:s}".format(fpath, time_stamp)
        os.rename(fpath, new_fpath)


def read_txt(file_txt):
    """ read in a text file as a string
    """
    with open(file_txt, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()

    return file_str
