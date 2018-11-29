""" Install automechanic
"""
from distutils.core import setup


setup(name="automechanic",
      version="0.1.1",
      packages=["automechanic",
                "automechanic.cli", "automechanic.task",
                "automechanic.mol", "automechanic.mol._irdkit",
                "automechanic.mol._ipyx2z", "automechanic.parse",
                "automechanic.parse", "automechanic.rere",
                # deprecated:
                "automechanic_old", "automechanic_old.ipybel",
                "automechanic_old.parse", "automechanic_old.parse.rere",
                "automechanic_old.ipyx2z", "from_qtc"],
      scripts=["automech",
               # deprecated:
               "automech_old"])
