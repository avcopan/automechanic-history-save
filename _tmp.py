""" testing script
"""
import automechanic.mol._irdkit as irdk
# from automechanic import mol

# mol.graph.conn.inchi(mol.inchi.connectivity_graph(
#     'InChI=1S/C2H2F2/c3-1-2-4/h1-2H'))

MLF = """
  automech  2D grid

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000 0.000 0.000 RAD=1 VAL=4
M  V30 2 C 1.000 0.000 0.000 RAD=1 VAL=4
M  V30 3 F 0.000 1.000 0.000 RAD=1 VAL=1
M  V30 4 CL 1.000 -1.000 0.000 RAD=1 VAL=1
M  V30 5 F 1.000 1.000 0.000 RAD=1 VAL=1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2 CFG=0
M  V30 2 1 1 3 CFG=0
M  V30 3 1 2 4 CFG=0
M  V30 4 1 2 5 CFG=0
M  V30 END BOND
M  V30 END CTAB
M  END
"""

RDM = irdk.from_molfile(MLF)
ICH = irdk.to_inchi(RDM)
print(ICH)
