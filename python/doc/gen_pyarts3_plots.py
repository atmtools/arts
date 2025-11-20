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
   :no-members:

Generic Plot Function
---------------------

.. autofunction:: plot

Plotting Functions by Type
---------------------------

""")

        for routine in routines:
            if routine in ["plot", "common", "matplotlib", "numpy"]:
                continue
            # Use the module name as a prefix in the title
            f.write(f"""
{routine}.plot
{'~' * (len(routine) + 5)}

.. autofunction:: pyarts3.plots.{routine}.plot

""")


generate(sys.argv[1])
