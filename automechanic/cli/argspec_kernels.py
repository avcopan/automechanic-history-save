""" library of command-line arguments for automech
"""

HELP = (
    'help',
    (
        ('action', 'store_true'),
        ('help', "show this help message and exit")
    )
)

SUBCMD = (
    'subcommand',
    (
        ('type', str),
    )
)
