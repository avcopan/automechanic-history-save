""" CLI commands
"""
import os
from . import arglib as al
from .arg import specifier_from_kernel as specifier
from .clihelp import call_subcommand
from .clihelp import call_task
from .. import task

RXN_CSV_DEF = 'reactions.csv'
SPC_CSV_DEF = 'species.csv'
RXN_CSV_CHAR = 'r'
SPC_CSV_CHAR = 's'

HOME_DIR = os.path.expanduser("~")
DB_NAME = 'automech_db'
DB_PREFIX_DEF = os.path.join(HOME_DIR, DB_NAME)
DB_PREFIX_CHAR = 'd'


def automech(argt):
    """ automech command
    """
    call_subcommand(
        argt,
        subcmds=(
            ('chemkin', chemkin),
            ('species', species),
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
    call_task(
        argt,
        task.chemkin.to_csv,
        specs=(
            specifier(
                al.MECHANISM_TXT, inp=True,
                extra_kwargs=(('nargs', '+'),),
                extra_helps=('followed by a CHEMKIN thermo file, if needed',),
                # value_map=(lambda x: tuple(map(os.path.abspath, x)))
                # ^ not necessary after deleting the PREFIX option
            ),
            specifier(
                al.REACTIONS_CSV, out=True, opt_char=RXN_CSV_CHAR.upper(),
                extra_kwargs=(('default', RXN_CSV_DEF),),
            ),
            specifier(
                al.SPECIES_CSV, out=True, opt_char=SPC_CSV_CHAR.upper(),
                extra_kwargs=(('default', SPC_CSV_DEF),),
            ),
        )
    )


def species(argt):
    """ species sub-command
    """
    call_subcommand(
        argt,
        subcmds=(
            ('expand_stereo', species_expand_stereo),
        )
    )


def species_expand_stereo(argt):
    """ fill in species guess geometries
    """
    call_task(
        argt,
        task.species.expand_stereo,
        specs=(
            specifier(
                al.GEOM_SPEC_KEY, inp=True,
                allowed_values=task.species.EXPAND_STEREO__SPC_ID_KEYS,
            ),
            specifier(
                al.SPECIES_CSV, inp=True,
            ),
            specifier(
                al.SPECIES_CSV, out=True, opt_char=SPC_CSV_CHAR.upper(),
                extra_kwargs=(('default', SPC_CSV_DEF),),
            ),
        )
    )
