""" I/O helpers
"""
from __future__ import unicode_literals
from builtins import open


def read_txt(file_txt):
    """ read in a text file as a string
    """
    with open(file_txt, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()

    return file_str
