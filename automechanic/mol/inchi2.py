""" functions operating on InChI strings
"""
from ._molfile2 import from_stereo_graph as _mlf_from_stereo_graph
from ._irdkit import from_molfile as _rdm_from_molfile
from ._irdkit import to_inchi as _rdm_to_inchi
from ._inchi_aux import numbering as _ich_aux_numbering


def from_stereo_graph_with_numbering(sgr):
    """ InChI string from a stereo graph
    """
    mlf = _mlf_from_stereo_graph(sgr)
    print(mlf)
    rdm = _rdm_from_molfile(mlf)
    ich, ich_aux = _rdm_to_inchi(rdm, with_aux_info=True)
    nums = _ich_aux_numbering(ich_aux)
    return ich, nums
