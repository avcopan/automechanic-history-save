""" functions for tracking the command-line argument vector
"""
import os
import sys
import logging
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from itertools import chain
from . import arglib as al
from .arg import specifier_key
from .arg import specifier_value_mapping
from .arg import interpret_specifier
from .arg import specifier_from_kernel
from .arg import set_specifier_keyword_value

LOG_NAME = specifier_from_kernel(al.LOG_NAME, opt_char='L', out=True)
LOG_LEVEL = specifier_from_kernel(al.LOG_LEVEL, opt_char='V')
PRINT_OUT = specifier_from_kernel(al.PRINT_OUT, opt_char='p')


def tracker(argv, pos):
    """ argument vector tracker
    """
    return (argv, pos)


def make_tracker(argv):
    """ make an argument vector tracker from sys.argv
    """
    argv[0] = os.path.basename(argv[0])
    return tracker(argv=argv, pos=0)


def increment_tracker(argt):
    """ advance the position of the argument vector tracker
    """
    argv, pos = argt
    pos += 1
    return tracker(argv=argv, pos=pos)


def call_name(argt):
    """ current call name
    """
    argv, pos = argt
    cmd_str = ' '.join(argv[:pos])
    return cmd_str


def command_line(argt):
    """ the full command line
    """
    argv, _ = argt
    cmd_str = ' '.join(argv[1:])
    return cmd_str


def has_arguments_left(argt):
    """ is this tracker at the end?
    """
    argv, pos = argt
    return pos < len(argv)


def current_position_argument(argt):
    """ current argument in the tracker
    """
    argv, pos = argt
    return argv[pos]


def last_subcommand(argt):
    """ get the last subcommand in the argument tracker
    """
    argv, pos = argt
    return argv[pos-1]


def log_file_name(argt):
    """ get a file name, based on the last subcommand of the argument tracker
    """
    base = last_subcommand(argt)
    fname = '{:s}.{:s}'.format(base, 'log')
    return fname


def current_argument_vector(argt):
    """ current command arguments
    """
    argv, pos = argt
    return argv[pos:]


def parse_current_argument(argt, spec):
    """ parse the argument at the current tracker position
    """
    cmd = call_name(argt)
    specs = [spec]
    par = _parser(cmd, specs)

    if (not has_arguments_left(argt) or
            current_position_argument(argt) in ('-h', '--help')):
        _quit_with_help_message(par)

    argv = [current_position_argument(argt)]
    val, = _parse_values(par, specs, argv)
    return val


def call_subcommand(argt, subcmds):
    """ parse argument to call appropriate subcommand function
    """
    argt = increment_tracker(argt)
    subcmd_keys, _ = zip(*subcmds)
    subcmd_key_sp = specifier_from_kernel(al.SUBCMD,
                                          allowed_values=subcmd_keys)
    subcmd_key_val = parse_current_argument(argt, subcmd_key_sp)
    subcmd_fnc = dict(subcmds)[subcmd_key_val]
    subcmd_fnc(argt)


def parse_arguments(argt, specs):
    """ parse the argument tracker according to specifications
    """
    cmd = call_name(argt)
    par = _parser(cmd, specs)

    argv = current_argument_vector(argt)

    if '-h' in argv or '--help' in argv:
        _quit_with_help_message(par)

    val_maps = tuple(map(specifier_value_mapping, specs))
    vals = _parse_values(par, specs, argv)
    vals = tuple(val_map(val) for val_map, val in zip(val_maps, vals))
    return vals


def call_task(argt, task, specs):
    """ parse remaining arguments to invoke a task
    """
    assert current_position_argument(argt) == task.__name__
    # ^ sanity check to keep the API and the CLI consistent with each other
    argt = increment_tracker(argt)
    fname = log_file_name(argt)
    log_name_spec = set_specifier_keyword_value(LOG_NAME, 'default', fname)

    extra_specs = (log_name_spec, LOG_LEVEL, PRINT_OUT)
    all_specs = tuple(chain(specs, extra_specs))
    all_vals = parse_arguments(argt, all_specs)

    nextra = len(extra_specs)
    vals, extra_vals = all_vals[:-nextra], all_vals[-nextra:]

    log_name, log_level, print_out = extra_vals
    logger = _logger(log_name, log_level, print_out)

    cmd_str = command_line(argt)
    logger.info("# command line: {:s}".format(cmd_str))

    task(*vals, logger=logger)


def _quit_with_help_message(par):
    par.print_help()
    par.exit()


def _parse_values(par, specs, argv):
    par_dct = vars(par.parse_args(argv))
    keys = tuple(map(specifier_key, specs))
    vals = tuple(map(par_dct.__getitem__, keys))
    return vals


def _parser(cmd, specs):
    par = ArgumentParser(
        prog=cmd,
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    for args, kwargs in map(interpret_specifier, specs):
        par.add_argument(*args, **kwargs)
    return par


def _logger(log_name, log_level, print_out):
    logger = logging.getLogger()
    logger.setLevel(log_level)

    fhandler = logging.FileHandler(log_name, mode='w')
    fhandler.setLevel(log_level)

    formatter = logging.Formatter('%(message)s')
    fhandler.setFormatter(formatter)

    logger.addHandler(fhandler)

    if print_out:
        shandler = logging.StreamHandler(sys.stdout)
        shandler.setLevel(log_level)
        shandler.setFormatter(formatter)
        logger.addHandler(shandler)

    return logger
