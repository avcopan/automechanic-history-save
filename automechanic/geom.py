""" cartesian geometry functions
"""
from .ipybel.geom import smiles
from .ipybel.geom import xyz_string
from .ipyx2z.geom import graph
from .ipyx2z.geom import resonance_graphs
from .ipyx2z.geom import radical_sites


def lewis_resonance_graphs(mgeo):
    mgrphs = resonance_graphs(mgeo)


def atoms(mgeo):
    """ the atoms in this cartesian geometry
    """
    atms, _ = zip(*mgeo)
    return atms


def formula(mgeo):
    """ the molecular formula of this cartesian geometry
    """
    atms = atoms(mgeo)
    return {atm: atms.count(atm) for atm in set(atms)}


__all__ = ['smiles', 'xyz_string', 'graph', 'resonance_graphs',
           'radical_sites']
