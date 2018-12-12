""" .xyz-based functions
"""
import re
from .parse_help import SPACE
from .parse_help import UNSIGNED_INTEGER
from .parse_help import STRING_START
from .parse_help import LINE_END
from .parse_help import LETTER
from .parse_help import FLOAT
from .parse_help import maybe
from .parse_help import one_or_more
from .parse_help import named_capture

SPACES = one_or_more(SPACE)


def number_of_atoms(dxyz):
    """ number of atoms from a .xyz string
    """
    natms_pattern = maybe(SPACES).join(
        [STRING_START, named_capture(UNSIGNED_INTEGER, 'natms'), LINE_END])
    match = re.search(natms_pattern, dxyz, re.MULTILINE)
    assert match
    gdct = match.groupdict()
    natms = int(gdct['natms'])
    return natms


def geometry(dxyz):
    """ geometry from a .xyz string
    """
    natms = number_of_atoms(dxyz)
    atomic_symbol = LETTER + maybe(LETTER)
    atom_pattern = SPACES.join(
        [named_capture(atomic_symbol, 'asymb'), named_capture(FLOAT, 'x'),
         named_capture(FLOAT, 'y'), named_capture(FLOAT, 'z')])
    line_pattern = atom_pattern + maybe(SPACES) + LINE_END

    mgeo = []
    for match in re.finditer(line_pattern, dxyz, re.MULTILINE):
        gdct = match.groupdict()
        asymb = gdct['asymb']
        xyz = tuple(map(float, [gdct['x'], gdct['y'], gdct['z']]))
        mgeo.append((asymb, xyz))

    if len(mgeo) != natms:
        raise ValueError("\nThis .xyz file is inconsistent: {:s}".format(dxyz))

    return tuple(mgeo)
