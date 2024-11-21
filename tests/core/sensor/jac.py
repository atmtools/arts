import pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
f = np.linspace(-5e6, 5e6, 1001) + line_f0

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequency(fmin=40e9, fmax=120e9, by_line=1)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_unit = "Tb"
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.ray_path_observer_agendaSet(option="Geometric")
ws.spectral_radiance_observer_agendaSet(option="Emission")

# %% Set up a sensor with Gaussian channel widths on individual frequency ranges

pos = [100e3, 0, 0]
los = [180.0, 0.0]

ws.measurement_sensorSimpleGaussian(
    frequency_grid=f, std=1e5, pos=pos, los=los, pol="Ih"
)
ws.measurement_sensorAddSimpleGaussian(
    frequency_grid=f, std=1e5, pos=pos, los=los, pol="Iv"
)
ws.measurement_sensorAddSimpleGaussian(
    frequency_grid=f, std=1e5, pos=pos, los=los, pol="RC"
)
ws.measurement_sensorAddSimpleGaussian(
    frequency_grid=f, std=1e5, pos=pos, los=los, pol="LC"
)

# %% Original calculations

ws.measurement_vectorFromSensor()
orig = ws.measurement_vector * 1.0

# %% Modify the sensor

DF1 = 1e6
DF2 = 2e6
DF3 = 3e6
DF4 = -4e6

ws.measurement_sensorSimpleGaussian(
    frequency_grid=f + DF1, std=1e5, pos=pos, los=los, pol="Ih"
)
ws.measurement_sensorAddSimpleGaussian(
    frequency_grid=f + DF2, std=1e5, pos=pos, los=los, pol="Iv"
)
ws.measurement_sensorAddSimpleGaussian(
    frequency_grid=f + DF3, std=1e5, pos=pos, los=los, pol="RC"
)
ws.measurement_sensorAddSimpleGaussian(
    frequency_grid=f + DF4, std=1e5, pos=pos, los=los, pol="LC"
)

ws.measurement_vectorFromSensor()
mod = ws.measurement_vector * 1.0

# %% Set up the retrieval

ws.RetrievalInit()
ws.RetrievalAddSensorFrequencyPolyFit(elem=0, d=1e3, matrix=np.diag(np.ones((1)) * 1e6))
ws.RetrievalFinalizeDiagonal()

ws.measurement_vector_fitted = []
ws.model_state_vector = []
ws.measurement_jacobian = [[]]
# %%
