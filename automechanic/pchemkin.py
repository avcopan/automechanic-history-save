""" functions for parsing CHEMKIN files
"""
import re
from re import escape
from .parse import maybe
from .parse import capture
from .parse import zero_or_more
from .parse import one_or_more
from .parse import one_of_these
from .parse import repeat_range
from .parse import named_capture
from .parse import group_dictionary
from .parse import ANY_CHAR
from .parse import WHITESPACE
from .parse import PLUS
from .parse import NON_NEWLINE
from .parse import NEWLINE
from .parse import FLOAT
from .parse import SIGN
from .parse import INTEGER
from .parse import STRING_START
from .parse import LINE_START


WHITESPACES = one_or_more(WHITESPACE, greedy=False)
ARROW = maybe(escape('<')) + escape('=') + maybe(escape('>'))
PADDED_PLUS = maybe(WHITESPACES) + PLUS + maybe(WHITESPACES)
PADDED_ARROW = maybe(WHITESPACES) + ARROW + maybe(WHITESPACES)
PADDED_EM = maybe(WHITESPACES) + 'M' + maybe(WHITESPACES)


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
    head = name + zero_or_more(WHITESPACE)
    foot = zero_or_more(WHITESPACE) + 'END'
    contents = one_or_more(ANY_CHAR, greedy=False)
    block_pattern = head + named_capture(contents, name='contents') + foot

    ret_str = None
    clean_mech_str = remove_comments(mech_str)
    gdct = group_dictionary(block_pattern, clean_mech_str)
    if gdct:
        ret_str = gdct['contents']
    return '\n'.join(ret_str.splitlines())


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


def reactions(mech_str, specs=None):
    """ all reactions in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: reactions
    :rtype: list of strings
    """
    specs = species(mech_str) if specs is None else specs
    reag_pattern = _en_reagents_pattern(specs)
    reac_pattern = _reaction_pattern(reag_pattern)
    reac_block_str = reactions_block(mech_str)
    return re.findall(reac_pattern, reac_block_str)


def therm_datas(mech_str, specs=None):
    """ all NASA polynomials in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: reactions
    :rtype: list of strings
    """
    polys = None

    specs = species(mech_str) if specs is None else specs
    ther_pattern = _therm_data_pattern(specs)
    ther_block_str = thermo_block(mech_str)
    if ther_block_str:
        polys = re.findall(ther_pattern, ther_block_str, re.MULTILINE)

    return polys


def split_reaction(reac, specs):
    """ split a CHEMKIN reaction into reactants and products

    :param reac: reaction string
    :type reac: str
    :param specs: species strings
    :type specs: list of strings

    :returns: reactants and products
    :rtype: (tuple of strings, tuple of strings)
    """
    nreag_pattern = capture(_en_reagents_pattern(specs))
    reactant_str, product_str = re.split(PADDED_ARROW, reac)
    reactant_ens = re.findall(nreag_pattern, reactant_str)
    product_ens = re.findall(nreag_pattern, product_str)
    reactants = sum((n * (spec,) for n, spec in
                     map(_split_en_reagents, reactant_ens)), ())
    products = sum((n * (spec,) for n, spec in
                    map(_split_en_reagents, product_ens)), ())
    return reactants, products


def split_therm_data(poly, specs):
    """ split a CHEMKIN NASA polynomial string into its constuent parts

    :param poly: NASA polynomial string
    :type poly: str
    :param specs: species strings
    :type specs: list of strings

    :returns: species, low T & high T coeffs, crossing T, low T & high T bounds
    :rtype: (str, list of floats, list of floats, float, float, float)
    """
    poly_data = None

    lines = poly.strip().splitlines()
    head_line = lines[0]
    coef_lines = '\n'.join(lines[1:])

    reag_pattern = named_capture(_reagent_pattern(specs), 'species')
    temps_pattern = WHITESPACES.join([named_capture(FLOAT, 'temp_lo'),
                                      named_capture(FLOAT, 'temp_hi'),
                                      named_capture(FLOAT, 'temp_cross')])

    pattern = WHITESPACES.join([reag_pattern, one_or_more(NON_NEWLINE),
                                temps_pattern])

    gdct = group_dictionary(pattern, head_line)

    exp = FLOAT + one_of_these(['E', 'e']) + maybe(SIGN) + INTEGER
    cfts = re.findall(exp, coef_lines)

    if gdct and len(cfts) in (14, 15):
        spec = gdct['species']
        temp_lo = float(gdct['temp_lo'])
        temp_hi = float(gdct['temp_hi'])
        temp_cross = float(gdct['temp_cross'])
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
    reagents = (reag_pattern + maybe(PADDED_PLUS + reag_pattern) + PADDED_PLUS
                + PADDED_EM)
    reac_pattern = reagents + PADDED_ARROW + reagents
    return reac_pattern


def _reaction_pattern_with_parentheses_em(reag_pattern):
    reagents = (reag_pattern + maybe(PADDED_PLUS + reag_pattern)
                + escape('(') + PADDED_PLUS + PADDED_EM + escape(')'))
    reac_pattern = reagents + PADDED_ARROW + reagents
    return reac_pattern


def _split_en_reagents(reag_str):
    count_pattern = STRING_START + maybe(capture(INTEGER))
    count_pattern_ = named_capture(count_pattern, name='count')
    gdct = group_dictionary(count_pattern_, reag_str)
    count = int(gdct['count']) if gdct['count'] else 1
    name = re.sub(count_pattern, '', reag_str)
    return count, name


def _en_reagents_pattern(specs):
    return maybe(INTEGER) + _reagent_pattern(specs)


def _reagent_pattern(specs):
    return one_of_these(_long_to_short(map(escape, specs)))


def _long_to_short(iterable):
    return list(reversed(sorted(iterable, key=len)))


def _therm_data_pattern(specs):
    reagent = _reagent_pattern(specs)
    maybe_spaces = maybe(WHITESPACES)
    exp = FLOAT + one_of_these(['E', 'e']) + maybe(SIGN) + INTEGER

    thermo_pattern = maybe_spaces.join(
        [LINE_START, reagent, WHITESPACE, one_or_more(NON_NEWLINE),
         FLOAT, WHITESPACE, FLOAT, WHITESPACE, FLOAT,
         zero_or_more(NON_NEWLINE),
         WHITESPACE, maybe(INTEGER), '1', NEWLINE,
         exp, exp, exp, exp, exp,
         WHITESPACE, '2', NEWLINE,
         exp, exp, exp, exp, exp,
         WHITESPACE, '3', NEWLINE,
         exp, exp, exp, exp, maybe(exp),
         WHITESPACE, '4', NEWLINE])
    return thermo_pattern
