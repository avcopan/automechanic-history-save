""" Geometry-based interface to PyBel
"""
import pybel


def smiles(geom):
    """ SMILES string associated with a cartesian geometry
    """
    pbmol = _pybel_molecule(geom)
    smi = pbmol.write('can').strip()
    return smi


def _pybel_molecule(geom):
    """ pybel.Molecule object of a cartesian geometry
    """
    dxyz = _xyz_string(geom)
    pbmol = pybel.readstring('xyz', dxyz)
    return pbmol


def _xyz_string(geom):
    """ .xyz format string of a cartesian geometry
    """
    natms = len(geom)
    dxyz = '{:d}\n\n'.format(natms)
    for asymb, xyz in geom:
        dxyz += '{:s} {:s} {:s} {:s}\n'.format(asymb, *map(repr, xyz))
    return dxyz
