""" command-line interface for automechanic
"""


def main(argv):
    """ automech main function
    """
    from .cmds import automech
    from .clihelp import make_tracker

    argt = make_tracker(argv)
    automech(argt)
