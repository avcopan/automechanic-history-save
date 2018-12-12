""" test the automechanic.fs module
"""
import tempfile
from automechanic import fs
from automechanic import fslib

C8H13O_MULT = 2
C8H13O_ICHS = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-/m0/s1',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4+/t8-/m0/s1',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3+,6-4-/t8-/m0/s1',
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3+,6-4+/t8-/m0/s1'
)


def test__branch__create():
    """ test fs.branch.create
    """
    fs_root_pth = tempfile.mkdtemp()

    with fs.enter(fs_root_pth):
        mult = C8H13O_MULT
        for ich in C8H13O_ICHS:
            sgms = fslib.species.branch_segments(ich, mult)
            fs.branch.create(sgms)
            fs.branch.validate(sgms)


if __name__ == '__main__':
    test__branch__create()
