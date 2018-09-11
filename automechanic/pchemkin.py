""" functions for parsing CHEMKIN files
"""
import re
from re import escape
from .parse import maybe
from .parse import zero_or_more
from .parse import one_or_more
from .parse import capture
from .parse import named_capture
from .parse import one_of_these
from .parse import repeat_range
from .parse import group_dictionary
from .parse import group_lists
from .parse import ANY_CHAR
from .parse import WHITESPACE
from .parse import PLUS
from .parse import NON_NEWLINE
from .parse import NEWLINE
from .parse import FLOAT
from .parse import SIGN
from .parse import INTEGER
from .parse import STRING_START


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


def thermo_lower_coefficients(mech_str):
    """ low-T NASA coefficients in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: dictionary, keyed by species name
    :rtype: dict
    """
    return _thermo(mech_str, 'cfts_lo')


def thermo_upper_coefficients(mech_str):
    """ high-T NASA coefficients in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: dictionary, keyed by species name
    :rtype: dict
    """
    return _thermo(mech_str, 'cfts_hi')


def thermo_common_temperatures(mech_str):
    """ low-T NASA coefficients in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: dictionary, keyed by species name
    :rtype: dict
    """
    return _thermo(mech_str, 'temp_md')


def reaction_strings(mech_str):
    """ all reaction strings in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: reactions
    :rtype: list
    """
    specs = species(mech_str)
    reag_pattern = _reagent_multiplet_pattern(specs)
    reac_pattern = _reaction_pattern_any(reag_pattern)
    reac_block_str = reactions_block(mech_str)
    return re.findall(reac_pattern, reac_block_str)


def reactions(mech_str):
    """ all reactions in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: reactions
    :rtype: list
    """
    specs = species(mech_str)
    reag_pattern = _reagent_multiplet_pattern(specs)
    reac_pattern = _reaction_pattern_any(reag_pattern)
    reacs = _reactions(mech_str, reac_pattern, specs)
    return reacs


def reactions_without_em(mech_str):
    """ pressure-inependent reactions in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: reactions
    :rtype: list
    """
    specs = species(mech_str)
    reag_pattern = _reagent_multiplet_pattern(specs)
    reac_pattern = _reaction_pattern_without_em(reag_pattern)
    reacs = _reactions(mech_str, reac_pattern, specs)
    return reacs


def reactions_with_em(mech_str):
    """ low-pressure limit reactions in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: reactions
    :rtype: list
    """
    specs = species(mech_str)
    reag_pattern = _reagent_multiplet_pattern(specs)
    reac_pattern = _reaction_pattern_with_em(reag_pattern)
    reacs = _reactions(mech_str, reac_pattern, specs)
    return reacs


def reactions_with_parentheses_em(mech_str):
    """ low-pressure limit reactions in a CHEMKIN file string

    :param mech_str: CHEMKIN file contents
    :type mech_str: str

    :returns: reactions
    :rtype: list
    """
    specs = species(mech_str)
    reag_pattern = _reagent_multiplet_pattern(specs)
    reac_pattern = _reaction_pattern_with_parentheses_em(reag_pattern)
    reacs = _reactions(mech_str, reac_pattern, specs)
    return reacs


def _reaction_pattern_any(reag_pattern):
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


def _reactions(mech_str, reac_pattern, specs):
    reac_block_str = reactions_block(mech_str)
    reac_strs = re.findall(reac_pattern, reac_block_str)
    reacs = tuple(_split_reaction(reac_str, specs) for reac_str in reac_strs)
    return reacs


def _split_reaction(reac_str, specs):
    nreag_pattern = capture(_reagent_multiplet_pattern(specs))
    rstr, pstr = re.split(PADDED_ARROW, reac_str)
    reactant_multiplets = tuple(map(_split_reagent_multiplet,
                                    re.findall(nreag_pattern, rstr)))
    product_multiplets = tuple(map(_split_reagent_multiplet,
                                   re.findall(nreag_pattern, pstr)))
    reactants = sum((cnt * (spec,) for cnt, spec in reactant_multiplets), ())
    products = sum((cnt * (spec,) for cnt, spec in product_multiplets), ())
    return reactants, products


def _split_reagent_multiplet(reag_str):
    count_pattern = STRING_START + maybe(capture(INTEGER))
    count_pattern_ = named_capture(count_pattern, name='count')
    gdct = group_dictionary(count_pattern_, reag_str)
    count = int(gdct['count']) if gdct['count'] else 1
    name = re.sub(count_pattern, '', reag_str)
    return count, name


def _reagent_multiplet_pattern(specs):
    return maybe(INTEGER) + _reagent_pattern(specs)


def _reagent_pattern(specs):
    return one_of_these(_long_to_short(map(escape, specs)))


def _long_to_short(iterable):
    return list(reversed(sorted(iterable, key=len)))


def _thermo(mech_str, key):
    """ gross
    """
    assert key in ('temp_lo', 'temp_hi', 'temp_md', 'cfts_lo', 'cfts_hi')

    specs = species(mech_str)
    reagent = _reagent_multiplet_pattern(specs)
    maybe_spaces = maybe(WHITESPACES)
    exp = FLOAT + one_of_these(['E', 'e']) + maybe(SIGN) + INTEGER

    thermo_pattern = maybe_spaces.join(
        [capture(reagent), WHITESPACE, one_or_more(NON_NEWLINE),
         capture(FLOAT), WHITESPACE, capture(FLOAT), WHITESPACE,
         capture(FLOAT), zero_or_more(NON_NEWLINE),
         WHITESPACE, '1', NEWLINE,
         capture(exp), capture(exp), capture(exp), capture(exp), capture(exp),
         WHITESPACE, '2', NEWLINE,
         capture(exp), capture(exp), capture(exp), capture(exp), capture(exp),
         WHITESPACE, '3', NEWLINE,
         capture(exp), capture(exp), capture(exp), capture(exp),
         WHITESPACE, '4', NEWLINE])

    ther_block_str = thermo_block(mech_str)
    grps = group_lists(thermo_pattern, ther_block_str)

    spec_lst = tuple(grp[0] for grp in grps)

    if key == 'temp_lo':
        data_lst = tuple(float(grp[1]) for grp in grps)
    elif key == 'temp_hi':
        data_lst = tuple(float(grp[2]) for grp in grps)
    elif key == 'temp_md':
        data_lst = tuple(float(grp[3]) for grp in grps)
    elif key == 'cfts_hi':
        data_lst = tuple(tuple(map(float, grp[4:11])) for grp in grps)
    elif key == 'cfts_lo':
        data_lst = tuple(tuple(map(float, grp[11:])) for grp in grps)

    return dict(zip(spec_lst, data_lst))
