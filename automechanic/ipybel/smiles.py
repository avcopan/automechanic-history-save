""" SMILES-based interface to PyBel
"""
import pybel


def canonical(smi):
    """ canonical smiles string
    """
    pbmol = _pybel_molecule(smi, geom=False)
    can_smi = pbmol.write('can').strip()
    return can_smi


def xyz_string(smi):
    """ XYZ string of the atoms in a SMILES string

    :param smi: SMILES string
    :type smi: str

    :returns: .xyz format string
    :rtype: string
    """
    pbmol = _pybel_molecule(smi)
    dxyz = pbmol.write('xyz').strip()
    return dxyz


def geometry(smi):
    """ cartesian geometry of the atoms in a SMILES string

    :param smi: SMILES string
    :type smi: str

    :returns: atomic symbols and their cartesian coordinates
    :rtype: list of pairs of strings and triples of floats
    """
    pbmol = _pybel_molecule(smi)
    anums = atomic_numbers(smi)
    asymbs = list(map(_atomic_symbol, anums))
    coords = tuple(tuple(atom.coords) for atom in pbmol.atoms)
    atms = tuple((asymb, xyz) for asymb, xyz in zip(asymbs, coords))
    return atms


def number_of_atoms(smi):
    """ number of atoms in a SMILES string

    :param smi: SMILES string
    :type smi: str

    :rtype: int
    """
    asymbs = atomic_symbols(smi)
    return len(asymbs)


def formula(smi):
    """ molecular formula of the atoms in a SMILES string

    :param smi: SMILES string
    :type smi: str

    :rtype: dict
    """
    asymbs = atomic_symbols(smi)
    return {asymb: asymbs.count(asymb) for asymb in set(asymbs)}


def atomic_symbols(smi):
    """ atomic symbols of the atoms in a SMILES string

    :param smi: SMILES string
    :type smi: str

    :rtype: list of strings
    """
    anums = atomic_numbers(smi)
    asymbs = list(map(_atomic_symbol, anums))
    return asymbs


def atomic_numbers(smi):
    """ atomic numbers of the atoms in a SMILES string

    :param smi: SMILES string
    :type smi: str

    :rtype: list of integers
    """
    pbmol = _pybel_molecule(smi, geom=False)
    anums = [atom.atomicnum for atom in pbmol.atoms]
    return anums


def _pybel_molecule(smi, geom=True):
    pbmol = pybel.readstring('smi', smi)
    pbmol.addh()
    if geom:
        pbmol.make3D()
    return pbmol


def _atomic_symbol(anum):
    """ convert atomic number to atomic symbol
    """
    anum2asymb = ['X', 'H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE',
                  'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA',
                  'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',
                  'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR',
                  'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN',
                  'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND',
                  'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB',
                  'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
                  'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH',
                  'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM',
                  'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS',
                  'RG', 'UUB', 'UUT', 'UUQ', 'UUP', 'UUH', 'UUS', 'UUO']
    return anum2asymb[anum]
