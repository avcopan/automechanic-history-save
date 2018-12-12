""" filesystem branch functions
"""
import os
from itertools import starmap as _starmap
from itertools import accumulate as _accumulate
from more_itertools import consume as _consume
from ..iohelp import write_yaml

INFO_FILE_NAME = 'info.yaml'


def create(sgms):
    """ creates a branch and returns the branch path
    """
    dir_names, inf_dcts = zip(*sgms)
    dir_pths = _dir_path_sequence(dir_names)
    _consume(map(_make_dir, dir_pths))
    _consume(_starmap(_write_info, zip(dir_pths, inf_dcts)))
    return dir_pths[-1]


def validate(sgms):
    """ assert that all directories and info files exist (for testing)
    """
    dir_names, inf_dcts = zip(*sgms)
    dir_pths = _dir_path_sequence(dir_names)
    _consume(map(_validate_dir, dir_pths))
    _consume(_starmap(_validate_info, zip(dir_pths, inf_dcts)))


def _info_file_name(dir_pth):
    return os.path.join(dir_pth, INFO_FILE_NAME)


def _dir_path_sequence(dir_names):
    return tuple(_accumulate(dir_names, os.path.join))


def _make_dir(dir_pth):
    if not os.path.exists(dir_pth):
        os.mkdir(dir_pth)


def _write_info(dir_pth, inf_dct):
    if inf_dct is not None:
        inf_pth = _info_file_name(dir_pth)
        write_yaml(inf_pth, inf_dct)


def _validate_dir(dir_pth):
    assert os.path.exists(dir_pth) and os.path.isdir(dir_pth)


def _validate_info(dir_pth, inf_dct):
    if inf_dct is not None:
        inf_pth = _info_file_name(dir_pth)
        assert os.path.exists(inf_pth) and os.path.isfile(inf_pth)
