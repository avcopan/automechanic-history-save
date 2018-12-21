""" stereo graph library; by analogy to the InChI stereo layer

vertices: atomic symbols, implicit hydrogen counts, stereo parity
    (('O', 1, None), ('C', 1, None), ('C', 1, False), ...)
edges: bond connectivity and stereo parity
    {{0, 1}: None, {1, 2}: True, ...}

stereo atom values:
    None = no stereo
    False = negative-parity stereo
    True = positive-parity stereo
stereo bond values:
    None = no stereo
    False = negative-parity stereo
    True = positive-parity stereo
"""
