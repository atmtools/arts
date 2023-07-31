# Import the module
import pyarts
import numpy as np  # This example uses numpy
import matplotlib.pyplot as plt


# Create a workspace
ws = pyarts.Workspace()

# The main agenda for computations of the radiative transfer
# follows a pure clears sky radiative transfer in this example.
# When the path tracing hits the surface, the surface emission
# is computed from the lowest atmospheric temperature,
# and when the path tracing hits the top of the atmosphere,
# the cosmic background is used as emission source.
ws.iy_main_agendaSet(option="Emission")
ws.iy_surface_agendaSet(option="UseSurfaceRtprop")
ws.surface_rtprop_agendaSet(option="Blackbody_SurfTFromt_field")
ws.iy_space_agendaSet(option="CosmicBackground")

# The path tracing is done step-by-step following a geometric path
ws.ppath_agendaSet(option="FollowSensorLosPath")
ws.ppath_step_agendaSet(option="GeometricPath")

# We might have to compute the Jacobian for the retrieval
# of relative humidity, so we need to set the agenda for
# that.
ws.water_p_eq_agendaSet(option="MK05")

# The geometry of our planet, and several other properties, are
# set to those of Earth.
ws.PlanetSet(option="Earth")

# Our output unit is Planc brightness temperature
ws.iy_unit = "PlanckBT"

# We do not care about polarization in this example
ws.stokes_dim = 1

# The atmosphere is assumed 1-dimensional
ws.atmosphere_dim = 1

# We have one sensor position in space, looking down, and one
# sensor position at the surface, looking up.
ws.sensor_pos = [[300e3], [0]]
ws.sensor_los = [[180], [0]]

# The dimensions of the problem are defined by a 1-dimensional pressure grid
ws.p_grid = np.logspace(5.01, -1)
ws.lat_grid = []
ws.lon_grid = []

# The surface is at 0-meters altitude
ws.z_surface = [[0.0]]

# Our sensor sees 10 frequency bins between 40 and 120 GHz
NF = 10
ws.f_grid = np.linspace(40e9, 120e9, NF)

# The atmosphere consists of water, oxygen and nitrogen.
# We set these to be computed using predefined absorption
# models.
ws.abs_speciesSet(species=["H2O-PWR98", "O2-PWR98", "N2-SelfContStandardType"])

# We have no line-by-line calculations, so we mark the LBL catalog as empty
ws.abs_lines_per_speciesSetEmpty()

# We now have all the information required to compute the absorption agenda.
ws.propmat_clearsky_agendaAuto()

# We need an atmosphere.  This is taken from the arts-xml-data in
# raw format before being turned into the atmospheric fields that
# ARTS uses internally.
ws.AtmRawRead(basename="planets/Earth/Fascod/tropical/tropical")
ws.AtmFieldsCalc()

# These calculations do no partial derivatives, so we can turn it off
ws.jacobianOff()

# There is no scattering in this example, so we can turn it off
ws.cloudboxOff()

# The concept of a sensor does not apply to this example, so we can turn it off
ws.sensorOff()

# We check the atmospheric geometry, the atmospheric fields, the cloud box,
# and the sensor.
ws.atmgeom_checkedCalc()
ws.atmfields_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()

# We perform the calculations
ws.yCalc()

# Create a simple plot to look at the simulations.  Try increasing NF
# above to see more details
plt.plot(ws.f_grid.value/1e9, ws.y.value.value.reshape(2, NF).T)
plt.xlabel("Frequency [GHz]")
plt.ylabel("Brightness temperature [K]")
plt.legend(["Looking down", "Looking up"])
plt.title("Low resolution O$_2$ millimeter absorption band")

# TESTING
# AS THIS FILE IS RUN TO TEST ARTS, WE NEED TO CHECK THAT
# THE CONTENT OF THE VARIABLES ARE GOOD.
assert np.allclose(
    ws.y.value,
    np.array([296.61207605, 292.04188873, 210.01486215, 271.44447592,
              293.06468686, 294.18513535, 293.97640914, 293.24271601,
              291.01085103, 252.58148054,  44.8272092 ,  88.45315869,
              292.69787616, 229.39778329, 104.89723532, 107.59101072,
              120.4062287 , 136.72275807, 160.08765173, 264.38696518])
    )
# END TESTING
