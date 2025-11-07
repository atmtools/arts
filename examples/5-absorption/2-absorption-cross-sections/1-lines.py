"""
This file will showcase how you can load line data from arts-cat-data into the
workspace and do the required setups to perform a simple forward calculations
using this data.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

# Initialize ARTS workspace
ws = pyarts.workspace.Workspace()

"""
The workspace contains an absorption species helper list
that allows it to understand what data and methods you
will want to call.

This example sets the absorption species to all oxygen
and all water isotopologues (in their line-by-line mode).
"""
ws.absorption_speciesSet(species=["O2-66", "H2O-161"])

"""
We now need to load the data.

ARTS comes with a lot of line-by-line data available
in its arts-cat-data repository.  The first line below
downloads that data to a local cache and allows you
to run this script.

You can also provide line data in other ways, e.g., by loading HITRAN
or some other format.  What is important is to populate the absorption
bands with appropriate data.
"""
pyarts.data.download()
ws.abs_bandsReadSpeciesSplitCatalog(basename="lines/")

"""
Compute absorption

Now we can use any of the methods available to compute line-by-line absorption.

Below we simply set a simple atmosphere and plot the resulting absorption.

Please see other examples for more details on how to use line-by-line
calculations in increasingly complex ways to solve your problem
"""

# Initialize an atmospheric point
atm = pyarts.arts.AtmPoint()
atm.temperature = 295  # At room temperature
atm.pressure = 1e5  # At 1 bar
atm["O2"] = 0.21  # At 21% atmospheric Oxygen
atm["H2O"] = 0.001  # At 0.1% atmospheric Water Vapor

# Set a frequency range and remove lines outside (to speed up calculations)
freqs = np.linspace(1e9, 1000e9, 1001)
ws.abs_bands.keep_frequencies(fmax=freqs[-1])

# Use available plotting routines to draw the results
fig, ax = pyarts.plot(
    ws.abs_bands,
    freqs=freqs,
    atm=atm,
    mode='important fill isotopes',
    min_pm=1e-7)
ax.legend(bbox_to_anchor=(1.05, 1))
ax.set_xlim(freqs[0], freqs[-1])
ax.set_ylim(ymin=1e-7, ymax=1)
ax.set_xticks(np.linspace(0, 1000e9, 11), np.linspace(0, 1000, 11))
ax.set_yscale('log')
ax.set_xlabel("Frequency [GHz]")
ax.set_ylabel("Absorption [1/m]")
ax.set_title("Absorption by species")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()
