"""

This file will showcase how you can load collision-induced absorption (CIA)
data from arts-cat-data into the workspace and do the required setups to
perform a simple forward calculations using this data

Note that this example presumes that you have set the environment variable
ARTS_DATA_PATH to contain a path to a local copy of both arts-cat-data and
arts-xml-data before you import pyarts.  Please check that this is the case
if the example does not work for you.  You can easily check if this path is
set by adding the following two lines at the top of this pyarts-controlfile:

```
import os
print(os.environ.get("ARTS_DATA_PATH"))
```

"""

import pyarts
import numpy as np
import matplotlib.pyplot as plt

# Initialize ARTS
ws = pyarts.workspace.Workspace()

"""

Set ws.abs_species to the species tags that you wish to use.  CIA tags are
structured so that the main and secondary species are separated by the word
"CIA"

This example sets up self CIA by molecular oxygen

CIA is not just self-induced but can occur between two different molecules.  In
this case you can change the tag to read something like "O2-CIA-N2", though
you must ensure that there also exists another N2 species in your ws.abs_species
as this is the only way to pass the volume mixing ratio of the atmosphere into
ARTS

"""
ws.absorption_speciesSet(species=["O2-CIA-O2"])

"""

Loads all CIA data from a given folder.  This command expects the file
"cia/O2-CIA-O2.xml" to be found either by relative local path or in the
user-defined search paths.

"""
ws.absorption_cia_dataReadSpeciesSplitCatalog(basename="cia/")


"""

This example does not deal with line-by-line absorption at all.  We must still
ensure that the line-by-line catalog has the correct size and that it has been
set so that our automatic agenda routine can do its work

"""
ws.absorption_bands = []

"""

You should generally always call this after you are done setting up your
ws.abs_species and ws.abs_lines_per_species.  It will deal with the internal
ARTS setup for you.  Note that the flag use_abs_lookup=1 can be passed to this
method call to set up the agenda for USING the the lookup-table.  Without the
flag, ARTS should be configured correctly to either COMPUTE the lookup-table
or to compute the absorption on-the-fly

In this case, it turns out that the temparature extrapolation is not enough
for a tropical atmosphere scenario below, so we extend it a small bit.  Play
with this "T_extrapolfac" value to see the relevant error message

"""
ws.propagation_matrix_agendaAuto(T_extrapolfac=1)

"""

Compute absorption

Now we can use the propagation_matrix_agenda to compute the absorption of O2-66.
We can also use this agenda in more complicated setups that might require
absorption calculations, but that is for other examples

To just execute the agenda we need to still define its both its inputs and the
inputs required to initialize the propagation matrix

"""

ws.jacobian_targets = pyarts.arts.JacobianTargets()
ws.frequency_grid = pyarts.arts.convert.wavelen2freq(np.linspace(6900e-9, 5900e-9, 1001))
ws.atmospheric_point.temperature = 295  # At room temperature
ws.atmospheric_point.pressure = 1e5  # At 1 bar
ws.atmospheric_point[ws.absorption_species[0]] = 0.21  # At 21% atmospheric Oxygen
ws.propagation_path_point # No particular POSLOS

# Call the agenda with inputs above
ws.propagation_matrix_agendaExecute()

# Plot the absorption of this example
plt.figure(1)
plt.clf()
plt.plot(
    1e9 * pyarts.arts.convert.freq2wavelen(ws.frequency_grid.value),
    ws.propagation_matrix[:, 0],
)
plt.xlabel("Wavelength [nm]")
plt.ylabel("Absorption [1/m]")
plt.title("O2-CIA-O2 absorption from examples/arts-cat-data/cia/cia.py")

"""
That's it!  You are done and have reached the end of this example.  Everything
below here is just to ensure that ARTS does not break in the future.  It can
be safely ignored

"""
# Save test results
# ws.propagation_matrix.savexml("cia_test_result.xml", type="ascii")

# test that we are still OK
propagation_matrix_agenda = \
    pyarts.arts.PropmatVector.fromxml("cia_test_result.xml")
assert np.allclose(
    propagation_matrix_agenda, ws.propagation_matrix
), "O2 Absorption has changed"
