"""

This file will showcase how you can load collision-induced absorption (CIA)
data from arts-cat-data into the workspace and do the required setups to
perform a simple forward calculations using this data

Note that this example presumes that you have set ARTS_DATA_PATH to contain
a local copy of both arts-cat-data and arts-xml-data.  Please check that this
is the case if the example does not work for you.  You can easily check if this
path is set by adding the following two lines at the top of this
pyarts-controlfile:

```
import os
print(os.environ.get("ARTS_DATA_PATH"))
```

"""

import pyarts
import numpy as np

# Initialize ARTS
ws = pyarts.workspace.Workspace()

"""

Set ws.abs_species to the species tags that you wish to use.  CIA tags are
structured so that the main and secondary species are surrounded by the word
"CIA".  These tags are followed by an extraction index, to select which CIA
data is used.

This sets up self CIA by molecular oxygen for the first "O2-CIA-O2-0" CIA
band.  To include all CIA band by O2, you have to write:
    "O2-CIA-O2-0,O2-CIA-O2-1,O2-CIA-O2-2,O2-CIA-O2-3,O2-CIA-O2-4,O2-CIA-O2-5,O2-CIA-O2-6,O2-CIA-O2-7"
    
CIA is not just self-induced but can occur between two different moleules.  In
this case you can change the tag to read something like "O2-CIA-N2-0", though
you must ensure that there also exist another N2 species in your ws.abs_species
as this is the only way to pass the volume mixing ratio of the atmosphere into
ARTS

"""
ws.abs_speciesSet(species=["O2-CIA-O2-0"])

"""

Loads all CIA data from a given folder.  This command expects the file
"cia/O2-CIA-O2.xml" to be found either by relative local path or in the
user-defined search paths.

"""
ws.abs_cia_dataReadSpeciesSplitCatalog(basename="cia/")


# We will not use line-absorption in this example, please see arts-cat-data/lines
# for an example of using line-by-line data in ARTS.  Still, the line data variable has
# to have the correct size, which the following call ensures
ws.abs_lines_per_speciesSetEmpty()

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
ws.propmat_clearsky_agendaAuto(T_extrapolfac=1)

# That's it! You now have working absorption calculations for your loaded lines.
# We still need to setup an atmosphere to perform our calculations in.
# This part is not covered in great details in this example

# setup some agendas and core settings
ws.iy_unit="PlanckBT"
ws.iy_main_agendaSet(option="Emission")
ws.ppath_agendaSet(option="FollowSensorLosPath")
ws.ppath_step_agendaSet(option="GeometricPath")
ws.water_p_eq_agendaSet()
ws.iy_space_agendaSet()
ws.iy_surface_agendaSet()
ws.surface_rtprop_agendaSet(option="Blackbody_SurfTFromt_surface")

# set the frequency range of these computations --- change to compute other frequencies
ws.f_grid = ws.abs_cia_data.value[0].data[0].grids[0]

# set the atmosphere and planet
ws.PlanetSet(option="Earth")
ws.AtmosphereSet1D()
ws.p_grid = np.logspace(5.05, 0, 40)
ws.lat_grid = []
ws.lon_grid = []
ws.AtmRawRead(basename="planets/Earth/Fascod/tropical/tropical")
ws.AtmFieldsCalc()
ws.z_surface = [[0]]
ws.t_surface = [[300]]

# Calculation setup:
ws.stokes_dim = 1
ws.jacobianOff()
ws.cloudboxOff()

# sensor setup
ws.sensor_pos = [[300e3]]
ws.sensor_los = [[180]]
ws.sensorOff()

# run checks to ensure quality of setup
ws.lbl_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()

# perform the forward calculations using standard methods
ws.yCalc()

# Save test results
# ws.y.value.savexml("cia_test_result.xml", type="zascii")

# test that we are still 
y = pyarts.arts.Vector()
y.readxml("cia_test_result.xml")
assert np.allclose(y, ws.y.value), "yCalc has changed output in test"
