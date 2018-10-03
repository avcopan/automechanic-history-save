""" functions for generating parsing command-line arguments
"""


def specifier(key, kwargs, opt_char=None, vals=()):
    """ CLI argument specifier
    """
    kwargs = freeze_specifier_kwargs(kwargs)
    return (key, kwargs, opt_char, vals)


def specifier_key(spec):
    """ specifier key
    """
    key, _, _, _ = spec
    return key


def specifier_kwargs(spec):
    """ specifier kwargs
    """
    _, kwargs, _, _ = spec
    return kwargs


def specifier_option_character(spec):
    """ character for optional argument
    """
    _, _, opt_char, _ = spec
    return opt_char


def specifier_allowed_values(spec):
    """ allowed values
    """
    _, _, _, vals = spec
    return vals


def freeze_specifier_kwargs(kwargs):
    """ make specifier kwargs immutable
    """
    return tuple(dict(kwargs).items())


def unfreeze_specifier_kwargs(kwargs):
    """ make specifier kwargs mutable
    """
    return dict(kwargs)


def set_specifier_kwarg(spec, spec_key, spec_val):
    """ set the value of a specifier kwarg
    """
    (key, kwargs, opt_char, vals) = spec
    kwargs = unfreeze_specifier_kwargs(kwargs)
    kwargs[spec_key] = spec_val
    return specifier(key, kwargs, opt_char, vals)


def specifier_from_kernel(kernel, opt_char=None, allowed_values=(),
                          extra_helps=()):
    """ specifier from specifier kernel
    """
    req_key_fmt = "<{:s}>".format

    key, kwargs = kernel

    if opt_char is None:
        key = req_key_fmt(key)

    spec = specifier(
        key=key,
        kwargs=kwargs,
        opt_char=opt_char,
        vals=allowed_values
    )
    help_msg = help_message(spec, helps=extra_helps, vals=allowed_values)
    spec = set_specifier_kwarg(spec, 'help', help_msg)
    return spec


def help_message(spec, helps, vals):
    """ add help messages to the specifier
    """
    kwargs = unfreeze_specifier_kwargs(specifier_kwargs(spec))

    if 'help' in kwargs:
        helps = (kwargs['help'],) + helps
    if vals:
        val_str = 'options: {:s}'.format(', '.join(vals))
        helps = tuple(helps) + (val_str,)

    help_msg = '; '.join(helps)
    return help_msg


def interpret_specifier(spec):
    """ map specifier to args and kwargs for argparse.add_argument
    """
    opt_char_fmt = "-{:s}".format
    opt_key_fmt = "--{:s}".format

    key = specifier_key(spec)
    kwargs = unfreeze_specifier_kwargs(specifier_kwargs(spec))
    opt_char = specifier_option_character(spec)

    if opt_char is None:
        args = ()
        kwargs['dest'] = key
    elif isinstance(opt_char, str) and len(opt_char) == 1:
        args = (opt_char_fmt(opt_char), opt_key_fmt(key))
    else:
        raise ValueError("Invalid option key.")

    return args, kwargs
