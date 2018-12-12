""" filesystem module

terminology:
    - filesystem: a collection of trees
    - tree: a collection of branches
    - branch: a collection of segments
    - segment: a pair `(dir_name, inf_dct)` where `dir_name` is a string and
      `inf_dct` is a dictionary, used to invoke the creation of a single
      directory (`dir_name/`) and an information file in the YAML format
      (`dir_name/inf.yaml`)
"""
from . import branch
from .context import enter

__all__ = ['branch', 'enter']
