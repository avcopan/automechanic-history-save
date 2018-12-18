""" stereo graph library; by analogy to the InChI stereo layer

vertices: atomic symbols, implicit hydrogen counts, stereo parity
    (('O', 1, None), ('C', 1, None), ('C', 1, None), ...)
edges: bond connectivity and stereo parity
    {{0, 1}: None, {1, 2}: True, ...}
"""
