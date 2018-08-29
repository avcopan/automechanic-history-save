""" Install automechanic
"""
from distutils.core import setup


setup(name="automechanic",
      version="0.1.0",
      packages=["automechanic", "automechanic.ipybel", "automechanic.ipyx2z",
                "from_qtc"],
      scripts=["automech", "automech_init", "automech_abstr_init",
               "automech_abstr_run"])
