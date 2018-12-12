""" I/O helpers
"""
from __future__ import unicode_literals
from builtins import open
import os
import time
import yaml


def timestamp_if_exists(file_pth):
    """ open a file, avoiding overwrites if requested
    """
    if os.path.isfile(file_pth):
        time_stamp = time.strftime("%Y%m%d-%H%M%S")
        new_file_pth = "{:s}_{:s}".format(file_pth, time_stamp)
        os.rename(file_pth, new_file_pth)


def read_string(file_pth):
    """ read in a file as a string
    """
    with open(file_pth, encoding='utf8', errors='ignore') as file_obj:
        string = file_obj.read()

    return string


# def write_string(file_pth, string):
#     """ write a string to a file
#     """
#     with open(file_pth, mode='w') as file_obj:
#         file_obj.write(string)


# def read_yaml(file_pth):
#     """ read in a yaml file as a dictionary
#     """
#     with open(file_pth) as file_obj:
#         dct = yaml.load(file_obj)
#         assert isinstance(dct, dict)
#
#     return dct


def write_yaml(file_pth, dct):
    """ write a dictionary to a yaml file
    """
    assert isinstance(dct, dict)
    with open(file_pth, mode='w') as file_obj:
        dct = yaml.dump(dct, file_obj, default_flow_style=False)

    return dct
