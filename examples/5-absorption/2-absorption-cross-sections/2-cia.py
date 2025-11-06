"""

This file will showcase how you can load collision-induced absorption (CIA)
data from arts-cat-data into the workspace and do the required setups to
perform a simple forward calculations using this data

Note that this example presumes that you have set the environment variable
ARTS_DATA_PATH to contain a path to a local copy of both arts-cat-data and
arts-xml-data before you import pyarts3.  Please check that this is the case
if the example does not work for you.  You can easily check if this path is
set by adding the following two lines at the top of this pyarts-controlfile:

```
import os
print(os.environ.get("ARTS_DATA_PATH"))
```

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

This example sets the absorption species to O2 collision-induced
absorption (CIA) with self and also with N2
"""
ws.absorption_speciesSet(species=["O2-CIA-O2", "O2-CIA-N2"])

"""
We now need to load the data.

ARTS comes with a lot of CIA data available
in its arts-cat-data repository.  The first line below
downloads that data to a local cache and allows you
to run this script.

You can also provide CIA data in other ways, e.g., by loading HITRAN
or some other format.  What is important is to populate the absorption
CIA data object with appropriate data.
"""
pyarts.data.download()
ws.absorption_cia_dataReadSpeciesSplitCatalog(basename="cia/")

"""
Compute absorption

Now we can use any of the methods available to compute CIA absorption.

Below we simply set a simple atmosphere and plot the resulting absorption.

Please see other examples for more details on how to use CIA
calculations in increasingly complex ways to solve your problem
"""

# Initialize an atmospheric point
atm = pyarts.arts.AtmPoint()
atm.temperature = 295  # At room temperature
atm.pressure = 1e5  # At 1 bar
atm["O2"] = 0.21  # At 21% atmospheric Oxygen
atm["N2"] = 0.78  # At 78% atmospheric Nitrogen

# Plot the absorption of this example
fig, ax = pyarts.plot(ws.absorption_cia_data, atm=atm)
for i, a in enumerate(ax.flatten()):
    a.set_xlabel("Wavelength [nm]")
    a.set_ylabel("Absorption [1/m]")
    f0 = np.inf
    f1 = -np.inf
    for x in ws.absorption_cia_data[i].data:
        f0 = min(f0, x.grids[0][0])
        f1 = max(f1, x.grids[0][-1])
    f = np.linspace(f0, f1, 7)
    a.set_xticks(f, (pyarts.arts.convert.freq2wavelen(f)*1e9).round(1))
    a.set_yscale('log')
    a.set_title(a.get_lines()[0].get_label())
fig.suptitle("CIA Absorption")

if "ARTS_HEADLESS" not in os.environ:
    plt.tight_layout()
    plt.show()
