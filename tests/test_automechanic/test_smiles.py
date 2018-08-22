""" test the automechanic.smiles module
"""
import numpy
from automechanic import smiles


def test__geometry():
    """ test smiles.geometry()
    """
    atms = smiles.geometry('C')
    asymbs, coords = zip(*atms)
    assert asymbs == ('C', 'H', 'H', 'H', 'H')
    coords_rmsd = numpy.sqrt(numpy.sum(numpy.var(coords, axis=0), axis=0))
    assert numpy.shape(coords) == (5, 3)
    assert numpy.isclose(coords_rmsd, 0.9768928266)
