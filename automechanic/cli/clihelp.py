""" functions for tracking the command-line argument vector
"""
import os
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from .argspec import specifier_key
from .argspec import interpret_specifier


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


def parse_arguments(argt, specs):
    """ parse the argument tracker according to specifications
    """
    cmd = call_name(argt)
    par = _parser(cmd, specs)

    argv = current_argument_vector(argt)

    if '-h' in argv or '--help' in argv:
        _quit_with_help_message(par)

    vals = _parse_values(par, specs, argv)
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


def _quit_with_help_message(par):
    par.print_help()
    par.exit()


def _parse_values(par, specs, argv):
    par_dct = vars(par.parse_args(argv))
    keys = tuple(map(specifier_key, specs))
    vals = tuple(map(par_dct.__getitem__, keys))
    return vals
