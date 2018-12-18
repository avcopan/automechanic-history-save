""" test the automechanic.mol module
"""
from automechanic import mol

AR_ICH = 'InChI=1S/Ar'

C2H2F2_SMI = 'F/C=C/F'
C2H2F2_SMI_NO_STEREO = 'FC=CF'
C2H2F2_ICH = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H/b2-1+'
C2H2F2_ICH_NO_STEREO = 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H'
C2H2F2_ICH_STEREO_UNKNOWN = 'InChI=1/C2H2F2/c3-1-2-4/h1-2H/b2-1?'
C2H2F2_GEO = (('F', (1.584822920001, -0.748486564300, -0.4271224303432)),
              ('C', (0.619219854789, 0.190165523008, -0.2716389279610)),
              ('C', (-0.635730620967, -0.1839138961594, -0.1803643636082)),
              ('F', (-1.602333181611, 0.736677675476, -0.02605091648865)),
              ('H', (0.916321356258, 1.229945559249, -0.2271265738271)),
              ('H', (-0.882300328471, -1.224388297273, -0.229635969682)))
C2H2F2_ICH_ORDER = (1, 2, 0, 3)

C8H13O_ICH = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-/m0/s1')
C8H13O_ICH_NO_ENANTIOMER = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/'
    'b5-3-,6-4-/t8-')
C8H13O_ICH_PARTIAL_STEREO = (
    'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3/b5-3-/t8-/m0/s1')
C8H13O_ICH_NO_STEREO = 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3'
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

C2H4CLF_SMI_NO_STEREO = 'C(Cl)(F)C'
C2H4CLF_ICH = 'InChI=1S/C2H4ClF/c1-2(3)4/h2H,1H3/t2-/m0/s1'
C2H4CLF_ICH_NO_STEREO = 'InChI=1S/C2H4ClF/c1-2(3)4/h2H,1H3'
C2H4CLF_ICH_INCOMPLETE_STEREO = 'InChI=1S/C2H4ClF/c1-2(3)4/h2H,1H3/t2-'
C2H4CLF_ICH_STEREO_UNKNOWN = 'InChI=1/C2H4ClF/c1-2(3)4/h2H,1H3/t2?'

# rdkit geometry failures:
RDKIT_FAIL_ICHS = (
    'InChI=1S/C6H10/c1-3-5-6-4-2/h3-6H,1-2H3',
    'InChI=1S/C5H10/c1-3-5-4-2/h3,5H,4H2,1-2H3',
    'InChI=1S/C4H7O2/c1-4(2)3-6-5/h1,3H2,2H3',
    'InChI=1S/C8H13O/c1-2-3-4-5-6-7-8-9/h2-3,6-7H,4-5,8H2,1H3')
# pybel geometry failures:
PYBEL_FAIL_ICHS = (
    'InChI=1S/C3H7O/c1-3(2)4/h3H,1-2H3',
    'InChI=1S/C3H7O2/c1-3(2)5-4/h3-4H,1H2,2H3',
    'InChI=1S/C3H7O4/c4-6-2-1-3-7-5/h4H,1-3H2',
    'InChI=1S/C3H6O3/c4-2-1-3-6-5/h2,5H,1,3H2',
    'InChI=1S/C3H7O4/c4-6-2-1-3-7-5/h1,4-5H,2-3H2',
    'InChI=1S/C3H6O3/c4-6-3-1-5-2-3/h3-4H,1-2H2',
    'InChI=1S/C3H6O/c1-2-4-3-1/h1-3H2',
    'InChI=1S/C3H2/c1-2-3-1/h1-2H',
    'InChI=1S/C4H9O/c1-3-4(2)5/h4H,3H2,1-2H3',
    'InChI=1S/C4H9O2/c1-3-4(2)6-5/h4H,3H2,1-2H3',
    'InChI=1S/C4H10O2/c1-3-4(2)6-5/h4-5H,3H2,1-2H3',
    'InChI=1S/C4H8O/c1-4-2-3-5-4/h4H,2-3H2,1H3',
    'InChI=1S/C4H8O/c1-2-4-5-3-1/h1-4H2',
    'InChI=1S/C4H8O/c1-3-4(2)5-3/h3-4H,1-2H3',
    'InChI=1S/C4H9O4/c1-4(8-6)2-3-7-5/h4-5H,2-3H2,1H3',
    'InChI=1S/C4H9O4/c5-7-3-1-2-4-8-6/h5H,1-4H2',
    'InChI=1S/C4H9O4/c1-3(7-5)4(2)8-6/h3-6H,1H2,2H3',
    'InChI=1S/C4H8O3/c1-3-4(7-5)2-6-3/h3-5H,2H2,1H3',
    'InChI=1S/C4H8O3/c5-7-2-1-4-3-6-4/h4-5H,1-3H2',
    'InChI=1S/C4H8O3/c5-7-4-1-2-6-3-4/h4-5H,1-3H2',
    'InChI=1S/C4H8O3/c1-3(7-5)4-2-6-4/h3-5H,2H2,1H3',
    'InChI=1S/C4H8O3/c1-3-4(7-3)2-6-5/h3-5H,2H2,1H3',
    'InChI=1S/C4H8O3/c1-4(7-6)2-3-5/h3-4,6H,2H2,1H3',
    'InChI=1S/C4H8O3/c5-3-1-2-4-7-6/h3,6H,1-2,4H2',
    'InChI=1S/C4H9O/c1-4(2)3-5/h4H,3H2,1-2H3',
    'InChI=1S/C4H9O2/c1-4(2,3)6-5/h1-3H3',
    'InChI=1S/C4H10O2/c1-4(2,3)6-5/h5H,1-3H3',
    'InChI=1S/C4H8O4/c1-4(8-6)2-7-3(4)5/h3,5-6H,2H2,1H3',
    'InChI=1S/C4H9O3/c1-2-4(5)3-7-6/h4,6H,2-3H2,1H3',
    'InChI=1S/C8H13O/c1-4-6-8(9)7(3)5-2/h4-8H,2H2,1,3H3',
    'InChI=1S/C4H7O3/c5-3-1-2-4-7-6/h1-2,6H,3-4H2',
    'InChI=1S/C4H7O2/c1-2-3-4-6-5/h3-4H,2H2,1H3',
    'InChI=1S/C4H9O3/c1-2-4(5)3-7-6/h4-5H,2-3H2,1H3',
    'InChI=1S/C4H6O/c1-2-4-3-5-4/h2,4H,1,3H2',
    'InChI=1S/C5H12O2/c1-3-4-5(2)7-6/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C5H12O2/c1-3-5(4-2)7-6/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C5H11O2/c1-3-4-5(2)7-6/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O2/c1-3-5(4-2)7-6/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O/c1-3-4-5(2)6/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O/c1-3-5(6)4-2/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O2/c1-3-5(4-2)7-6/h5-6H,1,3-4H2,2H3',
    'InChI=1S/C5H10O/c1-2-5-3-4-6-5/h5H,2-4H2,1H3',
    'InChI=1S/C5H10O/c1-5-3-2-4-6-5/h5H,2-4H2,1H3',
    'InChI=1S/C5H10O/c1-2-4-6-5-3-1/h1-5H2',
    'InChI=1S/C5H10O/c1-3-5-4(2)6-5/h4-5H,3H2,1-2H3',
    'InChI=1S/C5H10O/c1-4-3-5(2)6-4/h4-5H,3H2,1-2H3',
    'InChI=1S/C5H11O4/c1-2-5(9-7)3-4-8-6/h5-6H,2-4H2,1H3',
    'InChI=1S/C5H11O4/c1-5(9-7)3-2-4-8-6/h5-6H,2-4H2,1H3',
    'InChI=1S/C5H11O4/c6-8-4-2-1-3-5-9-7/h6H,1-5H2',
    'InChI=1S/C5H11O4/c1-4(8-6)3-5(2)9-7/h4-6H,3H2,1-2H3',
    'InChI=1S/C5H11O4/c1-2-5(9-7)3-4-8-6/h5,7H,2-4H2,1H3',
    'InChI=1S/C5H11O4/c1-3-5(9-7)4(2)8-6/h4-5,7H,3H2,1-2H3',
    'InChI=1S/C5H11O4/c1-5(9-7)3-2-4-8-6/h5-7H,1-4H2',
    'InChI=1S/C5H11O4/c1-3-5(9-7)4(2)8-6/h4-7H,1,3H2,2H3',
    'InChI=1S/C5H11O4/c1-2-5(9-7)3-4-8-6/h5-7H,1-4H2',
    'InChI=1S/C5H11O4/c1-2-3-5(9-7)4-8-6/h5-7H,1-4H2',
    'InChI=1S/C5H11O4/c6-8-4-2-1-3-5-9-7/h2,6-7H,1,3-5H2',
    'InChI=1S/C5H11O4/c1-2-3-5(9-7)4-8-6/h2,5-7H,3-4H2,1H3',
    'InChI=1S/C5H11O4/c1-5(9-7)3-2-4-8-6/h3,5-7H,2,4H2,1H3',
    'InChI=1S/C5H10O2/c1-3-4-5(2)7-6/h3,5-6H,1,4H2,2H3',
    'InChI=1S/C5H10O2/c1-3-4-5(2)7-6/h3-6H,1-2H3',
    'InChI=1S/C5H10O3/c1-2-4-5(8-6)3-7-4/h4-6H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c1-4(8-6)5-2-3-7-5/h4-6H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c1-4-2-5(8-6)3-7-4/h4-6H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c6-8-4-5-2-1-3-7-5/h5-6H,1-4H2',
    'InChI=1S/C5H10O3/c6-8-5-2-1-3-7-4-5/h5-6H,1-4H2',
    'InChI=1S/C5H10O3/c1-3-5(7-3)4(2)8-6/h3-6H,1-2H3',
    'InChI=1S/C5H10O3/c1-4-5(8-4)2-3-7-6/h4-6H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c1-2-5(8-7)3-4-6/h4-5,7H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c1-5(8-7)3-2-4-6/h4-5,7H,2-3H2,1H3',
    'InChI=1S/C5H10O3/c6-4-2-1-3-5-8-7/h4,7H,1-3,5H2',
    'InChI=1S/C5H10O3/c1-4(6)3-5(2)8-7/h5,7H,3H2,1-2H3',
    'InChI=1S/C5H10O3/c1-2-5(6)3-4-8-7/h7H,2-4H2,1H3',
    'InChI=1S/C5H10O3/c1-3-5(6)4(2)8-7/h4,7H,3H2,1-2H3',
    'InChI=1S/C5H9O2/c1-2-5(7)3-4-6/h4-5H,2-3H2,1H3',
    'InChI=1S/C5H9O2/c1-5(7)3-2-4-6/h4-5H,2-3H2,1H3',
    'InChI=1S/C5H9O2/c6-4-2-1-3-5-7/h4H,1-3,5H2',
    'InChI=1S/C5H9O2/c1-5(7)3-2-4-6/h2-4H2,1H3',
    'InChI=1S/C5H9O2/c1-2-5(7)3-4-6/h2-4H2,1H3',
    'InChI=1S/C5H9O2/c1-3-5(7)4(2)6/h4H,3H2,1-2H3',
    'InChI=1S/C5H8O2/c1-4(6)3-5(2)7/h3H2,1-2H3',
    'InChI=1S/C5H12O3/c1-4(6)3-5(2)8-7/h4-7H,3H2,1-2H3',
    'InChI=1S/C5H12O3/c1-2-5(8-7)3-4-6/h5-7H,2-4H2,1H3',
    'InChI=1S/C5H10O3/c1-4-3-5(2,6)8-7-4/h4,6H,3H2,1-2H3',
    'InChI=1S/C5H10O3/c1-2-5(6)3-4-7-8-5/h6H,2-4H2,1H3',
    'InChI=1S/C5H12O2/c1-4-5(2,3)7-6/h6H,4H2,1-3H3',
    'InChI=1S/C5H12O2/c1-4(2)5(3)7-6/h4-6H,1-3H3',
    'InChI=1S/C5H12O2/c1-5(2)3-4-7-6/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C5H11O2/c1-4-5(2,3)7-6/h4H2,1-3H3',
    'InChI=1S/C5H11O2/c1-4(2)5(3)7-6/h4-5H,1-3H3',
    'InChI=1S/C5H11O2/c1-5(2)3-4-7-6/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O/c1-4-5(2,3)6/h4H2,1-3H3',
    'InChI=1S/C5H11O/c1-4(2)5(3)6/h4-5H,1-3H3',
    'InChI=1S/C5H11O/c1-5(2)3-4-6/h5H,3-4H2,1-2H3',
    'InChI=1S/C5H11O2/c1-4(2)5(3)7-6/h4-6H,1H2,2-3H3',
    'InChI=1S/C5H10O/c1-5-2-3-6-4-5/h5H,2-4H2,1H3',
    'InChI=1S/C5H10O/c1-5(2)3-4-6-5/h3-4H2,1-2H3',
    'InChI=1S/C5H11O4/c1-4(8-6)5(2,3)9-7/h4,7H,1-3H3',
    'InChI=1S/C5H11O4/c1-5(2,9-7)3-4-8-6/h7H,3-4H2,1-2H3',
    'InChI=1S/C5H11O4/c1-5(4-9-7)2-3-8-6/h5-7H,1-4H2',
    'InChI=1S/C5H10O3/c6-8-2-1-5-3-7-4-5/h5-6H,1-4H2',
    'InChI=1S/C5H10O3/c1-4(8-6)5(2)3-7-5/h4,6H,3H2,1-2H3',
    'InChI=1S/C5H10O3/c1-5(4-7-5)2-3-8-6/h6H,2-4H2,1H3',
    'InChI=1S/C5H10O3/c6-8-4-5-1-2-7-3-5/h5-6H,1-4H2',
    'InChI=1S/C5H10O3/c1-5(2)4(8-5)3-7-6/h4,6H,3H2,1-2H3',
    'InChI=1S/C5H10O3/c1-5(4-8-6)2-3-7-5/h6H,2-4H2,1H3',
    'InChI=1S/C5H10O3/c1-5(2)4(8-6)3-7-5/h4,6H,3H2,1-2H3',
    'InChI=1S/C5H11O3/c1-4(8-7)5(2,3)6/h4,6H,1-3H3',
    'InChI=1S/C5H9O/c1-4-5(2,3)6/h4H,1H2,2-3H3',
    'InChI=1S/C5H9O3/c1-2-5(8-7)3-4-6/h4-5H,2-3H2,1H3',
    'InChI=1S/C5H9O3/c1-5(8-7)3-2-4-6/h4-5H,2-3H2,1H3',
    'InChI=1S/C5H9O3/c6-4-2-1-3-5-8-7/h4H,1-3,5H2',
    'InChI=1S/C5H10O/c1-3-5(6)4-2/h3-4H2,1-2H3',
    'InChI=1S/C5H9O/c1-3-5(6)4-2/h1,3-4H2,2H3',
    'InChI=1S/C5H9O3/c1-4(6)3-5(2)8-7/h5H,3H2,1-2H3',
    'InChI=1S/C5H9O3/c1-2-5(6)3-4-8-7/h2-4H2,1H3',
    'InChI=1S/C5H9O3/c1-3-5(6)4(2)8-7/h4H,3H2,1-2H3',
    'InChI=1S/C5H8O/c1-3-5(6)4-2/h3H,1,4H2,2H3',
    'InChI=1S/C5H11O3/c1-2-3-5(4-6)8-7/h5-6H,2-4H2,1H3',
    'InChI=1S/C5H11O3/c1-2-5(8-7)3-4-6/h5-6H,2-4H2,1H3',
    'InChI=1S/C5H11O3/c1-5(8-7)3-2-4-6/h5-6H,2-4H2,1H3',
    'InChI=1S/C5H11O3/c6-4-2-1-3-5-8-7/h6H,1-5H2',
    'InChI=1S/C6H14O2/c1-3-4-5-6(2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C6H14O2/c1-3-5-6(4-2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C6H13O2/c1-3-4-5-6(2)8-7/h6H,3-5H2,1-2H3',
    'InChI=1S/C6H13O2/c1-3-5-6(4-2)8-7/h6H,3-5H2,1-2H3',
    'InChI=1S/C6H13O/c1-3-4-5-6(2)7/h6H,3-5H2,1-2H3',
    'InChI=1S/C6H13O/c1-3-5-6(7)4-2/h6H,3-5H2,1-2H3',
    'InChI=1S/C6H13O2/c1-3-4-5-6(2)8-7/h6-7H,2-5H2,1H3',
    'InChI=1S/C6H12O/c1-2-3-6-4-5-7-6/h6H,2-5H2,1H3',
    'InChI=1S/C6H12O/c1-2-6-4-3-5-7-6/h6H,2-5H2,1H3',
    'InChI=1S/C6H12O/c1-6-4-2-3-5-7-6/h6H,2-5H2,1H3',
    'InChI=1S/C6H12O/c1-3-4-6-5(2)7-6/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C6H12O/c1-3-6-4-5(2)7-6/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C6H12O/c1-5-3-4-6(2)7-5/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C6H12O/c1-3-5-6(4-2)7-5/h5-6H,3-4H2,1-2H3',
    'InChI=1S/C6H13O4/c1-2-3-6(10-8)4-5-9-7/h6-7H,2-5H2,1H3',
    'InChI=1S/C6H13O4/c1-2-6(10-8)4-3-5-9-7/h6-7H,2-5H2,1H3',
    'InChI=1S/C6H13O4/c1-6(10-8)4-2-3-5-9-7/h6-7H,2-5H2,1H3',
    'InChI=1S/C6H13O4/c1-3-6(10-8)4-5(2)9-7/h5-7H,3-4H2,1-2H3',
    'InChI=1S/C6H13O4/c1-5(9-7)3-4-6(2)10-8/h5-7H,3-4H2,1-2H3',
    'InChI=1S/C6H13O4/c1-2-3-6(10-8)4-5-9-7/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O4/c1-3-4-6(10-8)5(2)9-7/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O4/c1-2-6(10-8)4-3-5-9-7/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H12O2/c1-3-4-5-6(2)8-7/h3,6-7H,1,4-5H2,2H3',
    'InChI=1S/C6H13O4/c1-5(9-7)3-4-6(2)10-8/h5-8H,1,3-4H2,2H3',
    'InChI=1S/C6H12O3/c1-2-3-6(9-8)4-5-7/h5-6,8H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-2-6(9-8)4-3-5-7/h5-6,8H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-6(9-8)4-2-3-5-7/h5-6,8H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-3-6(9-8)4-5(2)7/h6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H12O3/c1-5(7)3-4-6(2)9-8/h6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H12O3/c1-2-3-6(7)4-5-9-8/h8H,2-5H2,1H3',
    'InChI=1S/C6H12O3/c1-3-4-6(7)5(2)9-8/h5,8H,3-4H2,1-2H3',
    'InChI=1S/C6H12O3/c1-2-6(7)4-3-5-9-8/h8H,2-5H2,1H3',
    'InChI=1S/C6H12O3/c1-5(9-7)2-3-6-4-8-6/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-2-3-5-6(9-7)4-8-5/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-2-5(9-7)6-3-4-8-6/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-2-5-3-6(9-7)4-8-5/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-5(9-7)6-3-2-4-8-6/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-5-2-3-6(9-7)4-8-5/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c7-9-5-6-3-1-2-4-8-6/h6-7H,1-5H2',
    'InChI=1S/C6H12O3/c1-2-3-5-6(9-5)4-8-7/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-3-5(9-7)6-4(2)8-6/h4-7H,3H2,1-2H3',
    'InChI=1S/C6H12O3/c1-4(9-7)3-6-5(2)8-6/h4-7H,3H2,1-2H3',
    'InChI=1S/C6H12O3/c1-5-6(9-5)3-2-4-8-7/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-4-3-6(8-4)5(2)9-7/h4-7H,3H2,1-2H3',
    'InChI=1S/C6H12O3/c1-2-5-6(9-5)3-4-8-7/h5-7H,2-4H2,1H3',
    'InChI=1S/C6H12O3/c1-3-5-6(8-5)4(2)9-7/h4-7H,3H2,1-2H3',
    'InChI=1S/C4H7O/c1-4-2-3-5-4/h4H,1-3H2',
    'InChI=1S/C4H7O/c1-3-4(2)5-3/h3-4H,1H2,2H3',
    'InChI=1S/C5H9O/c1-2-5-3-4-6-5/h5H,1-4H2',
    'InChI=1S/C5H9O/c1-5-3-2-4-6-5/h5H,1-4H2',
    'InChI=1S/C5H9O/c1-4-3-5(2)6-4/h4-5H,1,3H2,2H3',
    'InChI=1S/C6H11O2/c1-2-3-6(8)4-5-7/h5-6H,2-4H2,1H3',
    'InChI=1S/C6H11O2/c1-2-6(8)4-3-5-7/h5-6H,2-4H2,1H3',
    'InChI=1S/C6H11O2/c1-6(8)4-2-3-5-7/h5-6H,2-4H2,1H3',
    'InChI=1S/C6H11O2/c1-3-6(8)4-5(2)7/h6H,3-4H2,1-2H3',
    'InChI=1S/C6H11O2/c1-6(8)4-2-3-5-7/h2-5H2,1H3',
    'InChI=1S/C6H11O2/c1-2-3-6(8)4-5-7/h2-5H2,1H3',
    'InChI=1S/C6H11O2/c1-3-4-6(8)5(2)7/h5H,3-4H2,1-2H3',
    'InChI=1S/C6H11O2/c1-3-5(7)6(8)4-2/h5H,3-4H2,1-2H3',
    'InChI=1S/C6H11O2/c1-3-6(8)4-5(2)7/h5H,3-4H2,1-2H3',
    'InChI=1S/C6H11O2/c1-2-6(8)4-3-5-7/h2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-2-3-6(7)4-5-9-8/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-2-6(7)4-3-5-9-8/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-6(7)4-2-3-5-9-8/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-2-3-4-6(5-7)9-8/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-3-4-6(7)5(2)9-8/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-3-6(7)4-5(2)9-8/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-5(7)3-4-6(2)9-8/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-6(9-8)4-2-3-5-7/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-2-3-6(9-8)4-5-7/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-3-4-6(9-8)5(2)7/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-3-5(7)6(4-2)9-8/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-3-6(9-8)4-5(2)7/h5-6,8H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-2-6(9-8)4-3-5-7/h6,8H,2-5H2,1H3',
    'InChI=1S/C6H12O/c1-3-4-5-6(2)7/h3-5H2,1-2H3',
    'InChI=1S/C6H13O3/c1-2-3-4-6(7)5-9-8/h6-7H,2-5H2,1H3',
    'InChI=1S/C6H13O3/c1-3-4-6(9-8)5(2)7/h5-7H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-3-4-6(7)5(2)9-8/h5-7H,3-4H2,1-2H3',
    'InChI=1S/C6H13O3/c1-3-5(7)6(4-2)9-8/h5-7H,3-4H2,1-2H3',
    'InChI=1S/C6H11O/c1-3-4-5-6(2)7/h3,6H,1,4-5H2,2H3',
    'InChI=1S/C7H16O2/c1-3-4-5-6-7(2)9-8/h7-8H,3-6H2,1-2H3',
    'InChI=1S/C7H16O2/c1-3-5-6-7(4-2)9-8/h7-8H,3-6H2,1-2H3',
    'InChI=1S/C7H16O2/c1-3-5-7(9-8)6-4-2/h7-8H,3-6H2,1-2H3',
    'InChI=1S/C7H15O2/c1-3-4-5-6-7(2)9-8/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O2/c1-3-5-6-7(4-2)9-8/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O2/c1-3-5-7(9-8)6-4-2/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O/c1-3-4-5-6-7(2)8/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O/c1-3-5-6-7(8)4-2/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O/c1-3-5-7(8)6-4-2/h7H,3-6H2,1-2H3',
    'InChI=1S/C7H15O2/c1-3-5-7(9-8)6-4-2/h7-8H,1,3-6H2,2H3',
    'InChI=1S/C7H14O/c1-2-3-4-7-5-6-8-7/h7H,2-6H2,1H3',
    'InChI=1S/C7H14O/c1-2-4-7-5-3-6-8-7/h7H,2-6H2,1H3',
    'InChI=1S/C7H14O/c1-2-7-5-3-4-6-8-7/h7H,2-6H2,1H3',
    'InChI=1S/C7H14O/c1-3-4-5-7-6(2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H14O/c1-3-4-7-5-6(2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H14O/c1-3-7-5-4-6(2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H14O/c1-6-4-3-5-7(2)8-6/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H14O/c1-3-5-7-6(4-2)8-7/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H14O/c1-3-6-5-7(4-2)8-6/h6-7H,3-5H2,1-2H3',
    'InChI=1S/C7H15O4/c1-2-3-4-7(11-9)5-6-10-8/h7-8H,2-6H2,1H3',
    'InChI=1S/C7H15O4/c1-2-4-7(11-9)5-3-6-10-8/h7-8H,2-6H2,1H3')


def test__geom__inchi_with_order():
    """ test mol.geom.inchi_with_order
    """
    assert (mol.geom.inchi_with_order(C2H2F2_GEO)
            == (C2H2F2_ICH, C2H2F2_ICH_ORDER))


def test__geom__inchi():
    """ test mol.geom.inchi
    """
    assert mol.geom.inchi(C2H2F2_GEO) == C2H2F2_ICH


def test__smiles__inchi():
    """ test mol.smiles.inchi
    """
    assert mol.smiles.inchi(C2H2F2_SMI) == C2H2F2_ICH


def test__molfile__inchi():
    """ test mol.molfile.inchi
    """
    mlf = ('\n  automech  2D grid\n\n'
           '  0  0  0  0  0  0  0  0  0  0999 V3000\n'
           'M  V30 BEGIN CTAB\n'
           'M  V30 COUNTS 2 1 0 0 0\n'
           'M  V30 BEGIN ATOM\n'
           'M  V30 1 C 0.000 0.000 0.000 RAD=2 VAL=3 CFG=0\n'
           'M  V30 2 F 0.000 0.000 0.000 RAD=1 VAL=1 CFG=0\n'
           'M  V30 END ATOM\n'
           'M  V30 BEGIN BOND\n'
           'M  V30 1 1 1 2\n'
           'M  V30 END BOND\n'
           'M  V30 END CTAB\n'
           'M  END\n')
    assert mol.molfile.inchi(mlf) == 'InChI=1S/CH2F/c1-2/h1H2'


def test__inchi__smiles():
    """ test mol.inchi.smiles
    """
    assert mol.smiles.inchi(mol.inchi.smiles(AR_ICH)) == AR_ICH
    assert (mol.smiles.inchi(mol.inchi.smiles(C8H13O_ICH_NO_STEREO))
            == C8H13O_ICH_NO_STEREO)
    assert (tuple(map(mol.smiles.inchi, map(mol.inchi.smiles, C8H13O_ICHS)))
            == C8H13O_ICHS)


def test__inchi__recalculate():
    """ test mol.inchi.recalculate
    """
    assert mol.inchi.recalculate(C2H2F2_ICH_NO_STEREO) == C2H2F2_ICH_NO_STEREO
    assert (mol.inchi.recalculate(C2H2F2_ICH_NO_STEREO, force_stereo=True)
            == C2H2F2_ICH_STEREO_UNKNOWN)


def test__inchi__is_closed():
    """ test mol.inchi.is_closed
    """
    assert mol.inchi.is_closed(C8H13O_ICH) is True
    assert mol.inchi.is_closed(C8H13O_ICH_PARTIAL_STEREO) is False
    assert mol.inchi.is_closed(C8H13O_ICH_NO_STEREO) is True
    assert mol.inchi.is_closed(C8H13O_ICH_NO_ENANTIOMER) is False


def test__inchi__prefix():
    """ test mol.inchi.prefix
    """
    assert mol.inchi.prefix(C2H2F2_ICH) == 'InChI=1S'
    assert mol.inchi.prefix(C2H2F2_ICH_NO_STEREO) == 'InChI=1S'
    assert mol.inchi.prefix(C2H2F2_ICH_STEREO_UNKNOWN) == 'InChI=1'


def test__inchi__version():
    """ test mol.inchi.version
    """
    assert mol.inchi.version(C2H2F2_ICH) == '1S'
    assert mol.inchi.version(C2H2F2_ICH_NO_STEREO) == '1S'
    assert mol.inchi.version(C2H2F2_ICH_STEREO_UNKNOWN) == '1'


def test__inchi__formula_layer():
    """ test mol.inchi.formula_layer
    """
    assert mol.inchi.formula_layer(C2H2F2_ICH) == 'C2H2F2'
    assert (mol.inchi.formula_layer('InChI=1S/2C2H5.Zn/c2*1-2;/h2*1H2,2H3;')
            == '2C2H5.Zn')


def test__inchi__sublayer():
    """ test mol.inchi.sublayer
    """
    assert mol.inchi.sublayer(C2H2F2_ICH, 'c') == '3-1-2-4'
    assert mol.inchi.sublayer(C2H2F2_ICH, 'h') == '1-2H'
    assert mol.inchi.sublayer(C2H2F2_ICH, 'b') == '2-1+'
    assert mol.inchi.sublayer(C2H2F2_ICH_STEREO_UNKNOWN, 'c') == '3-1-2-4'
    assert mol.inchi.sublayer(C2H2F2_ICH_STEREO_UNKNOWN, 'h') == '1-2H'
    assert mol.inchi.sublayer(C2H2F2_ICH_STEREO_UNKNOWN, 'b') == '2-1?'
    assert mol.inchi.sublayer(C2H2F2_ICH_NO_STEREO, 'c') == '3-1-2-4'
    assert mol.inchi.sublayer(C2H2F2_ICH_NO_STEREO, 'h') == '1-2H'
    assert mol.inchi.sublayer(C2H2F2_ICH_NO_STEREO, 'b') is None


def test__inchi__with_sublayers():
    """ test mol.inchi.with_sublayers
    """
    assert mol.inchi.with_sublayers(AR_ICH, ('c', 'h')) == AR_ICH
    assert (mol.inchi.with_sublayers(C2H2F2_ICH, ('c', 'h'))
            == 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H')
    assert (mol.inchi.with_sublayers(C2H2F2_ICH_NO_STEREO, ('c', 'h'))
            == 'InChI=1S/C2H2F2/c3-1-2-4/h1-2H')
    assert (mol.inchi.with_sublayers(C2H2F2_ICH_STEREO_UNKNOWN, ('c', 'h'))
            == 'InChI=1/C2H2F2/c3-1-2-4/h1-2H')
    assert (mol.inchi.with_sublayers(C8H13O_ICH_NO_STEREO, ('c', 'h'))
            == 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3')
    assert (mol.inchi.with_sublayers(C8H13O_ICH, ('c', 'h'))
            == 'InChI=1S/C8H13O/c1-3-5-7-8(9)6-4-2/h3-6,8H,7H2,1-2H3')


def test__inchi__compatible_stereoisomers():
    """ test mol.inchi.compatible_stereoisomers
    """
    assert (mol.inchi.compatible_stereoisomers(C8H13O_ICH_NO_STEREO)
            == C8H13O_ICHS)
    assert mol.inchi.compatible_stereoisomers(C8H13O_ICH) == (C8H13O_ICH,)


def test__inchi__inchi_key():
    """ test mol.inchi.inchi_key
    """
    assert mol.inchi.inchi_key(C2H2F2_ICH) == 'WFLOTYSKFUPZQB-OWOJBTEDSA-N'
    assert (mol.inchi.inchi_key(C2H2F2_ICH_NO_STEREO)
            == 'WFLOTYSKFUPZQB-UHFFFAOYSA-N')


def test__inchi__geometry():
    """ test mol.inchi.geometry
    """
    # make sure these run
    mol.inchi.geometry(C2H2F2_ICH)
    for ich in C8H13O_ICHS:
        mol.inchi.geometry(ich)
    for ich in RDKIT_FAIL_ICHS:
        mol.inchi.geometry(ich)
    for ich in PYBEL_FAIL_ICHS:
        mol.inchi.geometry(ich)


def test__inchi__connectivity_graph():
    """ test mol.inchi.connectivity_graph
    """
    assert (mol.inchi.connectivity_graph(C2H2F2_ICH)
            == ((('C', 1), ('C', 1), ('F', 0), ('F', 0)),
                {frozenset({0, 1}): None, frozenset({0, 2}): None,
                 frozenset({1, 3}): None}))


def test__inchi__stereo_graph():
    """ test mol.inchi.stereo_graph
    """
    print(mol.inchi.stereo_graph(C8H13O_ICH))


def test__inchi__has_unknown_stereo_elements():
    """ test mol.inchi.has_unknown_stereo_elements
    """
    assert (mol.inchi.has_unknown_stereo_elements(C8H13O_ICH)
            is False)
    assert (mol.inchi.has_unknown_stereo_elements(C8H13O_ICH_PARTIAL_STEREO)
            is True)
    assert (mol.inchi.has_unknown_stereo_elements(C8H13O_ICH_NO_STEREO)
            is True)
    assert (mol.inchi.has_unknown_stereo_elements(C8H13O_ICH_NO_ENANTIOMER)
            is False)


def test__inchi_key__is_standard_neutral():
    """ test mol.inchi_key.is_standard_neutral()
    """
    assert (mol.inchi_key.is_standard_neutral(
        mol.inchi.inchi_key(C2H2F2_ICH)) is True)
    assert (mol.inchi_key.is_standard_neutral(
        mol.inchi.inchi_key(C2H2F2_ICH_NO_STEREO)) is True)
    assert (mol.inchi_key.is_standard_neutral(
        mol.inchi.inchi_key(C2H2F2_ICH_STEREO_UNKNOWN)) is False)


def test__inchi_key__first_hash():
    """ test mol.inchi_key.first_hash()
    """
    assert (mol.inchi_key.first_hash(
        mol.inchi.inchi_key(C2H2F2_ICH)) == 'WFLOTYSKFUPZQB')
    assert (mol.inchi_key.first_hash(
        mol.inchi.inchi_key(C2H2F2_ICH_NO_STEREO)) == 'WFLOTYSKFUPZQB')
    assert (mol.inchi_key.first_hash(
        mol.inchi.inchi_key(C2H2F2_ICH_STEREO_UNKNOWN)) == 'WFLOTYSKFUPZQB')


def test__inchi_key__second_hash():
    """ test mol.inchi_key.second_hash()
    """
    assert (mol.inchi_key.second_hash(
        mol.inchi.inchi_key(C2H2F2_ICH)) == 'OWOJBTED')
    assert (mol.inchi_key.second_hash(
        mol.inchi.inchi_key(C2H2F2_ICH_NO_STEREO)) == 'UHFFFAOY')
    assert (mol.inchi_key.second_hash(
        mol.inchi.inchi_key(C2H2F2_ICH_STEREO_UNKNOWN)) == 'HXYFBOIP')


if __name__ == '__main__':
    # test__smiles__inchi()
    # test__molfile__inchi()
    # test__inchi__smiles()
    # # test__inchi__recalculate()
    # test__inchi__is_closed()
    # # test__inchi__prefix()
    # # test__inchi__formula_layer()
    # # test__inchi__sublayer()
    # test__inchi__with_sublayers()
    # # test__inchi__inchi_key()
    # # test__inchi__geometry()
    # # test__inchi__connectivity_graph()
    # test__inchi__has_unknown_stereo_elements()
    # test__inchi_key__is_standard_neutral()
    # test__inchi_key__first_hash()
    # test__inchi_key__second_hash()
    # # test__inchi__compatible_stereoisomers()
    # # test__geom__inchi_with_order()
    # # test__geom__atoms()
    # # test__geom__bonds()
    # # test__geom__graph()
    # test__inchi__connectivity_graph()
    # test__inchi__stereo_graph()
    test__inchi__prefix()
    test__inchi__version()
    test__inchi__formula_layer()
