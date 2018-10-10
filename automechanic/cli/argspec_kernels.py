""" library of command-line arguments for automech
"""

SUBCMD = (
    'subcommand',
    (
        ('type', str),
    )
)

MECHANISM_TXT = (
    'mechanism_txt',
    (
        ('type', str),
        ('help', "CHEMKIN mechanism file")
    )
)

REACTIONS_CSV = (
    'reactions_csv',
    (
        ('type', str),
        ('help', "CSV with reaction data")
    )
)

SPECIES_CSV = (
    'species_csv',
    (
        ('type', str),
        ('help', "CSV with species data")
    )
)

PREFIX = (
    'prefix',
    (
        ('type', str),
        ('default', '.'),
        ('help', "prefix for output")
    )
)

LOG_NAME = (
    'log_name',
    (
        ('type', str),
        ('help', "log file for automech output")
    )
)

LOG_LEVEL = (
    'log_level',
    (
        ('type', int),
        ('default', 20),
        ('help', "log verbosity: 20=INFO, 10=DEBUG")
    )
)

PRINT_OUT = (
    'print_out',
    (
        ('action', 'store_true'),
        ('help', "print log output to screen")
    )
)
