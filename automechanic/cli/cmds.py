""" CLI commands
"""
import os
from .clihelp import call_subcommand
from .clihelp import call_routine
from .argspec import specifier_from_kernel as specifier
from .argspec_kernels import MECHANISM_TXT
from .argspec_kernels import REACTIONS_CSV
from .argspec_kernels import SPECIES_CSV
from .. import routines

REACTIONS_CSV_NAME = 'reactions.csv'
SPECIES_CSV_NAME = 'species.csv'


def automech(argt):
    """ automech command
    """
    call_subcommand(
        argt,
        subcmds=(
            ('chemkin', chemkin),
        )
    )


def chemkin(argt):
    """ chemkin sub-command
    """
    call_subcommand(
        argt,
        subcmds=(
            ('to_csv', chemkin_to_csv),
        )
    )


def chemkin_to_csv(argt):
    """ parse CHEMKIN mechanism
    """
    call_routine(
        argt,
        routines.chemkin.to_csv,
        specs=(
            specifier(
                MECHANISM_TXT, inp=True,
                extra_kwargs=(('nargs', '+'),),
                extra_helps=('followed by a CHEMKIN thermo file, if needed',),
                value_map=(lambda x: tuple(map(os.path.abspath, x)))
            ),
            specifier(
                REACTIONS_CSV, out=True, opt_char='R',
                extra_kwargs=(('default', REACTIONS_CSV_NAME),),
            ),
            specifier(
                SPECIES_CSV, out=True, opt_char='S',
                extra_kwargs=(('default', SPECIES_CSV_NAME),),
            ),
        )
    )
