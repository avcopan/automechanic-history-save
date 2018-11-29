""" test the automechanic.mol module
"""
import pytest
from automechanic import mol

FCCF_SMI = 'F/C=C/F'
FCCF_SMI_NO_STEREO = 'FC=CF'
FCCF_ICH = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+'
FCCF_ICH_NO_STEREO = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H'
FCCF_ICH_STEREO_UNKNOWN = 'InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1?'
FCCF_GEO = (('F', (1.584822920001, -0.748486564300, -0.4271224303432)),
            ('C', (0.619219854789, 0.190165523008, -0.2716389279610)),
            ('C', (-0.635730620967, -0.1839138961594, -0.1803643636082)),
            ('F', (-1.602333181611, 0.736677675476, -0.02605091648865)),
            ('H', (0.916321356258, 1.229945559249, -0.2271265738271)),
            ('H', (-0.882300328471, -1.224388297273, -0.229635969682)))
FCCF_ICH_ORDER = (1, 2, 0, 3)


def test__ge__inchi_with_order():
    """ test mol.ge.inchi_with_order
    """
    assert (mol.ge.inchi_with_order(FCCF_GEO, strict=True)
            == (FCCF_ICH, FCCF_ICH_ORDER))


def test__ge__inchi():
    """ test mol.ge.inchi
    """
    assert mol.ge.inchi(FCCF_GEO, strict=True) == FCCF_ICH


def test__ic__inchi():
    """ test mol.ic.inchi
    """
    with pytest.raises(Exception):
        mol.ic.inchi(FCCF_ICH_NO_STEREO)

    assert mol.ic.inchi(FCCF_ICH_NO_STEREO, strict=False) == FCCF_ICH_NO_STEREO
    assert (mol.ic.inchi(FCCF_ICH_NO_STEREO, force_stereo=True, strict=True)
            == FCCF_ICH_STEREO_UNKNOWN)


def test__ic__inchi_key():
    """ test mol.ic.inchi_key
    """
    with pytest.raises(Exception):
        mol.ic.inchi_key(FCCF_ICH_NO_STEREO)

    assert mol.ic.inchi_key(FCCF_ICH) == 'WFLOTYSKFUPZQB-OWOJBTEDSA-N'
    assert (mol.ic.inchi_key(FCCF_ICH_NO_STEREO, strict=False)
            == 'WFLOTYSKFUPZQB-UHFFFAOYSA-N')


def test__ic__geometry():
    """ test mol.ic.geometry
    """
    # for now, just make sure these run
    mol.ic.geometry(FCCF_ICH)
    mol.ic.geometry(FCCF_ICH_NO_STEREO, strict=False)
    # and this one doesn't
    with pytest.raises(Exception):
        mol.ic.geometry(FCCF_ICH_NO_STEREO, strict=True)
        # ^ checks that the InChI-Key of the geometry matches the original
        # InChI string


def test__sm__inchi():
    """ test mol.sm.inchi
    """
    assert mol.sm.inchi(FCCF_SMI) == FCCF_ICH
    assert mol.sm.inchi(FCCF_SMI_NO_STEREO, strict=False) == FCCF_ICH_NO_STEREO
    assert (mol.sm.inchi(FCCF_SMI_NO_STEREO, force_stereo=True)
            == FCCF_ICH_STEREO_UNKNOWN)
    with pytest.raises(Exception):
        mol.sm.inchi(FCCF_SMI_NO_STEREO)


if __name__ == '__main__':
    test__ge__inchi_with_order()
    test__ic__inchi()
    test__ic__inchi_key()
    test__ic__geometry()
    test__sm__inchi()
