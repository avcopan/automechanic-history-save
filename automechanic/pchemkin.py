""" functions for parsing CHEMKIN files
"""
import re
from re import escape
from itertools import chain
from more_itertools import windowed
from .parse import maybe
from .parse import capture
from .parse import zero_or_more
from .parse import one_or_more
from .parse import one_of_these
from .parse import repeat_range
from .parse import named_capture
from .parse import group_lists
from .parse import group_dictionary
from .parse import ANY_CHAR
from .parse import SPACE
from .parse import PLUS
from .parse import FLOAT
from .parse import SIGN
from .parse import INTEGER
from .parse import STRING_START
from .parse import STRING_END


SPACES = one_or_more(SPACE, greedy=False)
ARROW = maybe(escape('<')) + escape('=') + maybe(escape('>'))
PADDED_PLUS = maybe(SPACES) + PLUS + maybe(SPACES)
PADDED_ARROW = maybe(SPACES) + ARROW + maybe(SPACES)
PADDED_EM = maybe(SPACES) + 'M' + maybe(SPACES)
PLUS_EM = PADDED_PLUS + PADDED_EM
PAREN_PLUS_EM = escape('(') + PLUS_EM + escape(')')


def remove_comments(mech_str):
    """ remove comments from a CHEMKIN filestring

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: commentr-free CHEMKIN filestring
    :rtype: str
    """
    comment = '!.*'
    clean_mech_str = re.sub(comment, '', mech_str)
    return clean_mech_str


def block(mech_str, name):
    """ block from a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str
    :param name: the block header ('ELEM', 'SPEC', 'THER', 'REAC')
    :type name: str

    :returns: block string
    :rtype: str
    """
    head = name + zero_or_more(SPACE)
    foot = zero_or_more(SPACE) + 'END'
    contents = one_or_more(ANY_CHAR, greedy=False)
    block_pattern = head + named_capture(contents, name='contents') + foot

    ret_str = None
    clean_mech_str = remove_comments(mech_str)
    gdct = group_dictionary(block_pattern, clean_mech_str)
    if gdct:
        ret_str = gdct['contents']
        ret_str = '\n'.join(ret_str.splitlines())
    return ret_str


def elements_block(mech_str):
    """ elements block from a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: elements block string
    :rtype: str
    """
    name = one_of_these(['ELEMENTS', 'ELEM'])
    return block(mech_str, name=name)


def species_block(mech_str):
    """ species block from a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: species block string
    :rtype: str
    """
    name = one_of_these(['SPECIES', 'SPEC'])
    return block(mech_str, name=name)


def thermo_block(mech_str):
    """ thermo block from a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: thermo block string
    :rtype: str
    """
    name = one_of_these(['THERMO ALL', 'THERMO', 'THER'])
    return block(mech_str, name=name)


def reactions_block(mech_str):
    """ reactions block from a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: reactions block string
    :rtype: str
    """
    name = one_of_these(['REACTIONS', 'REAC'])
    return block(mech_str, name=name)


def species(mech_str):
    """ species in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: species
    :rtype: list
    """
    spec_block_str = species_block(mech_str)
    specs = tuple(spec_block_str.split())
    return specs


def reactions(mech_str):
    """ all reactions in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: reactions
    :rtype: list of strings
    """
    specs = species(mech_str)
    reag_pattern = _en_reagents_pattern(specs)
    reac_pattern = _reaction_pattern(reag_pattern)
    reac_block_str = reactions_block(mech_str)
    reacs = [reac.strip() for reac in re.findall(reac_pattern, reac_block_str)]
    return reacs


def therm_data_strings(mech_str):
    """ all NASA polynomials in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: reactions
    :rtype: list of strings
    """
    polys = None

    ther_block_str = thermo_block(mech_str)
    if ther_block_str:
        ther_lines = tuple(map(str.strip, ther_block_str.splitlines()))
        polys = []
        for stanza_lines in windowed(ther_lines, 4):
            match = all(lin.endswith('{:d}'.format(num))
                        for lin, num in zip(stanza_lines, range(1, 5)))
            if match:
                stanza = '\n'.join(stanza_lines)
                polys.append(stanza)

    return polys


def split_reaction(reac):
    """ split a CHEMKIN reaction into reactants and products

    :param reac: reaction string
    :type reac: str

    :returns: reactants and products
    :rtype: (tuple of strings, tuple of strings)
    """
    em_pattern = one_of_these([PAREN_PLUS_EM + STRING_END,
                               PLUS_EM + STRING_END])

    reactant_str, product_str = re.split(PADDED_ARROW, reac)
    reactant_str = re.sub(em_pattern, '', reactant_str)
    product_str = re.sub(em_pattern, '', product_str)
    en_reactants = tuple(map(_expand_en_reagents,
                             map(str.strip,
                                 re.split(PADDED_PLUS, reactant_str))))
    en_products = tuple(map(_expand_en_reagents,
                            map(str.strip,
                                re.split(PADDED_PLUS, product_str))))
    reactants = tuple(chain(*en_reactants))
    products = tuple(chain(*en_products))
    return reactants, products


def split_therm_data(poly):
    """ split a CHEMKIN NASA polynomial string into its constuent parts

    :param poly: NASA polynomial string
    :type poly: str

    :returns: species, low T & high T coeffs, crossing T, low T & high T bounds
    :rtype: (str, list of floats, list of floats, float, float, float)
    """
    poly_data = None

    lines = poly.strip().splitlines()
    head_line = lines[0]
    coef_lines = '\n'.join(lines[1:])

    spec = head_line.split()[0]

    temps_pattern = SPACES.join([named_capture(FLOAT, 'temp_lo'),
                                 named_capture(FLOAT, 'temp_hi'),
                                 named_capture(FLOAT, 'temp_cross')])
    temp_groups = group_lists(temps_pattern, head_line)

    exp = FLOAT + one_of_these(['E', 'e']) + maybe(SIGN) + INTEGER
    cfts = re.findall(exp, coef_lines)

    if temp_groups and len(cfts) in (14, 15):
        temps = map(float, temp_groups[-1])
        temp_lo, temp_hi, temp_cross = temps
        cfts_hi = tuple(map(float, cfts[:7]))
        cfts_lo = tuple(map(float, cfts[7:14]))
        poly_data = (spec, cfts_lo, cfts_hi, temp_cross, temp_lo, temp_hi)

    return poly_data


def _reaction_pattern(reag_pattern):
    reac_pattern = one_of_these([
        _reaction_pattern_with_parentheses_em(reag_pattern),
        _reaction_pattern_with_em(reag_pattern),
        _reaction_pattern_without_em(reag_pattern)])
    return reac_pattern


def _reaction_pattern_without_em(reag_pattern):
    reagents = reag_pattern + repeat_range(PADDED_PLUS + reag_pattern, 0, 2)
    reac_pattern = reagents + PADDED_ARROW + reagents
    return reac_pattern


def _reaction_pattern_with_em(reag_pattern):
    reagents = (reag_pattern + maybe(PADDED_PLUS + reag_pattern) + PLUS_EM)
    reac_pattern = reagents + PADDED_ARROW + reagents
    return reac_pattern


def _reaction_pattern_with_parentheses_em(reag_pattern):
    reagents = (reag_pattern + maybe(PADDED_PLUS + reag_pattern) +
                PAREN_PLUS_EM)
    reac_pattern = reagents + PADDED_ARROW + reagents
    return reac_pattern


def _expand_en_reagents(reag_str):
    count_pattern = STRING_START + maybe(capture(INTEGER))
    count_pattern_ = named_capture(count_pattern, name='count')
    gdct = group_dictionary(count_pattern_, reag_str)
    count = int(gdct['count']) if gdct['count'] else 1
    name = re.sub(count_pattern, '', reag_str)
    return (name,) * count


def _en_reagents_pattern(specs):
    return maybe(INTEGER) + _reagent_pattern(specs)


def _reagent_pattern(specs):
    return one_of_these(_long_to_short(map(escape, specs)))


def _long_to_short(iterable):
    return list(reversed(sorted(iterable, key=len)))