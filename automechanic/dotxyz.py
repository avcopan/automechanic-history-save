""" .xyz-based functions
"""
import re
from .parse import SPACE
from .parse import INTEGER
from .parse import STRING_START
from .parse import LINE_START
from .parse import LINE_END
from .parse import LETTER
from .parse import FLOAT
from .parse import maybe
from .parse import one_or_more
from .parse import named_capture

SPACES = one_or_more(SPACE)


def number_of_atoms(dxyz):
    """ number of atoms from a .xyz string
    """
    natms_pattern = maybe(SPACES).join(
        [STRING_START, named_capture(INTEGER, 'natms'), LINE_END])
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
    line_pattern = LINE_START + atom_pattern + LINE_END

    mgeo = []
    for match in re.finditer(line_pattern, dxyz, re.MULTILINE):
        gdct = match.groupdict()
        asymb = gdct['asymb']
        xyz = tuple(map(float, [gdct['x'], gdct['y'], gdct['z']]))
        mgeo.append((asymb, xyz))

    if len(mgeo) != natms:
        raise ValueError("\nThis .xyz file is inconsistent: {:s}".format(dxyz))

    return tuple(mgeo)
