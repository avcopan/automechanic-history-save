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
FILESYSTEM_NAME = 'automech_fs'
FILESYSTEM_PREFIX_DEF = os.path.join(HOME_DIR, FILESYSTEM_NAME)
FILESYSTEM_PREFIX_CHAR = 'f'

STEREO_HANDLING_CHAR = 't'


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
            ('to_csv', chemkin__to_csv),
        )
    )


def chemkin__to_csv(argt):
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
            ('to_inchi', species__to_inchi),
            ('filesystem', species__filesystem),
        )
    )


def species__to_inchi(argt):
    """ expand species into their possible stereoisomers
    """
    call_task(
        argt,
        task.species.to_inchi,
        specs=(
            specifier(
                al.SPECIES_ID,
                allowed_values=task.species.TO_INCHI__SPC_ID_VALS,
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


def species__filesystem(argt):
    """ fill in species guess geometries
    """
    call_task(
        argt,
        task.species.filesystem,
        specs=(
            specifier(
                al.SPECIES_CSV, inp=True,
            ),
            specifier(
                al.SPECIES_CSV, out=True, opt_char=SPC_CSV_CHAR.upper(),
                extra_kwargs=(('default', SPC_CSV_DEF),),
            ),
            specifier(
                al.STEREO_HANDLING, opt_char=STEREO_HANDLING_CHAR,
                allowed_values=task.species.FILESYSTEM__STEREO_HANDLING_VALS,
                extra_kwargs=(
                    ('default', task.species.STEREO_HANDLING_DEF),),
            ),
            specifier(
                al.FILESYSTEM_PREFIX, out=True,
                opt_char=FILESYSTEM_PREFIX_CHAR.upper(),
                extra_kwargs=(('default', FILESYSTEM_PREFIX_DEF),),
            ),
        )
    )
