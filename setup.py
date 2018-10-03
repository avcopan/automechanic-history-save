""" Install automechanic
"""
from distutils.core import setup


setup(name="automechanic",
      version="0.2a1",
      packages=["automechanic", "automechanic.ipybel", "automechanic.ipyx2z",
                "automechanic.cli", "from_qtc"],
      scripts=["automech", "automech2"])
