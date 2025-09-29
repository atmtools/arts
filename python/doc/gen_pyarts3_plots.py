import sys
import pyarts3


def generate(path):
    routines = [r for r in dir(pyarts3.plots) if not r.startswith("_")]
    routines.sort()

    with open(f"{path}/pyarts3.plots.rst", "w") as f:
        f.write("""pyarts3.plots
=============

.. currentmodule:: pyarts3.plots

.. automodule:: pyarts3.plots
   :members:
   :imported-members:
   :undoc-members:
   :special-members: __init__, __call__

""")

        for routine in routines:
            f.write(f"""
.. rubric:: {routine}

.. automodule:: pyarts3.plots.{routine}
   :members:
   :imported-members:
   :undoc-members:
   :special-members: __init__, __call__
""")


generate(sys.argv[1])
