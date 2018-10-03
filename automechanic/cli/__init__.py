""" command-line interface for automechanic
"""
from .argspec_kernels import SUBCMD
from .argspec import specifier_from_kernel as specifier
from .clihelp import make_tracker
from .clihelp import parse_arguments
from .clihelp import parse_current_argument
from .clihelp import increment_tracker


def main(argv):
    """ automech main function
    """
    argt = make_tracker(argv)
    automech(argt)


def automech(argt):
    """ automech command
    """
    argt = increment_tracker(argt)
    subcmd_keys, subcmd_fncs = zip(
        ('chemkin', chemkin),
    )
    subcmd_key_sp = specifier(SUBCMD, allowed_values=subcmd_keys)
    subcmd_key_val = parse_current_argument(argt, subcmd_key_sp)
    subcmd_dct = dict(zip(subcmd_keys, subcmd_fncs))
    subcmd_fnc = subcmd_dct[subcmd_key_val]
    subcmd_fnc(argt)


def chemkin(argt):
    """ chemkin sub-command
    """
    argt = increment_tracker(argt)
    subcmd_keys, subcmd_fncs = zip(
        ('to_csv', chemkin_to_csv),
    )
    subcmd_key_sp = specifier(SUBCMD, allowed_values=subcmd_keys)
    subcmd_key_val = parse_current_argument(argt, subcmd_key_sp)
    subcmd_dct = dict(zip(subcmd_keys, subcmd_fncs))
    subcmd_fnc = subcmd_dct[subcmd_key_val]
    subcmd_fnc(argt)


def chemkin_to_csv(argt):
    """ parse CHEMKIN mechanism
    """
    argt = increment_tracker(argt)
    print "HI from chemkin_to_csv"
