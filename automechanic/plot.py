""" plotting functions
"""
from matplotlib import pyplot


def add_line(diag, x_vals, y_vals, dotted=False):
    """ add line plot to diagram
    """
    _, axs = diag
    lnsty = '--' if dotted else '-'
    axs.plot(x_vals, y_vals, linestyle=lnsty)


def add_axis_labels(diag, x_label, y_label):
    """ add axis labels to the diagram
    """
    _, axs = diag
    axs.set_xlabel(x_label)
    axs.set_ylabel(y_label)
    return diag


def add_text(diag, text_str, left=True, top=True, margin=0.05):
    """ add text to diagram
    """
    xpos = margin if left else 1. - margin
    xaln = 'left' if left else 'right'
    ypos = margin if not top else 1. - margin
    yaln = 'bottom' if not top else 'top'

    _, axs = diag
    axs.annotate(text_str, xy=(xpos, ypos), xycoords='axes fraction',
                 ha=xaln, va=yaln)


def write_diagram(diag, fname, close=False):
    """ write diagram to a file
    """
    fig, _ = diag
    fig.savefig(fname)
    if close:
        pyplot.close(fig)


def empty_diagram():
    """ canvas = (figure, axes)
    """
    return pyplot.subplots()
