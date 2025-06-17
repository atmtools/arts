import pyarts
import numpy as np
import matplotlib.pyplot as plt

LIMIT = 50
noise = 0.5
NF = 1001
std = 1
pol1 = "RC"
pol2 = "Ih"
RTOL = 5e-3

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
f = np.linspace(-5e6, 5e6, NF) + line_f0

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_unit = "Tb"
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.ray_path_observer_agendaSetGeometric()
ws.spectral_radiance_observer_agendaSet(option="Emission")

# %% Set up a sensor with Gaussian channel widths on individual frequency ranges

pos = [100e3, 0, 0]
los = [180.0, 0.0]

ws.measurement_sensorSimpleGaussian(
    frequency_grid=f, std=std, pos=pos, los=los, pol=pol1
)
ws.measurement_sensorAddSimpleGaussian(
    frequency_grid=f, std=std, pos=pos, los=los, pol=pol2
)

# %% Original calculations

ws.measurement_vectorFromSensor()
orig = ws.measurement_vector * 1.0

ws.measurement_sensorSimpleGaussian(
    frequency_grid=f, std=std, pos=pos, los=los, pol=pol1
)
ws.measurement_sensorAddSimpleGaussian(
    frequency_grid=f, std=std, pos=pos, los=los, pol=pol2
)

ws.measurement_vectorFromSensor()
mod = ws.measurement_vector * 1.0


# %% Set up the retrieval

ws.RetrievalInit()
ws.RetrievalAddErrorPolyFit(
    sensor_elem=0, t=f, matrix=np.diag(np.ones((1)) * 25), polyorder=0
)
ws.RetrievalAddErrorPolyFit(
    sensor_elem=1, t=f - min(f), matrix=np.diag([25, 1e-12]), polyorder=1
)
ws.RetrievalFinalizeDiagonal()

# %% Perform OEM retrieval of frequency grid

fail = True
for i in range(LIMIT):
    ws.measurement_vector_fitted = []
    ws.model_state_vector = []
    ws.measurement_jacobian = [[]]

    ws.measurement_sensorSimpleGaussian(
        frequency_grid=f, std=std, pos=pos, los=los, pol=pol1
    )
    ws.measurement_sensorAddSimpleGaussian(
        frequency_grid=f, std=std, pos=pos, los=los, pol=pol2
    )
    ws.measurement_vectorFromSensor()
    ws.measurement_vector += np.random.normal(0, noise, 2 * NF)
    ws.measurement_vector[:NF] += 5 + f * 0
    ws.measurement_vector[NF:] += 15 + (f - min(f)) * 1e-6

    ws.measurement_vector_fitted = []
    ws.model_state_vector = []
    ws.measurement_jacobian = [[]]

    ws.model_state_vector_aprioriFromData()

    ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)

    ws.OEM(method="gn")

    print(f"got:      {ws.model_state_vector:B,}")
    print(f"expected: [5, 15, 1e-6]")
    if np.allclose(ws.model_state_vector / [5, 15, 1e-6], 1, rtol=RTOL):
        print(f"Within {RTOL}%.  Success!")
        fail = False
        break
    print(f"RelDiff not less than {RTOL}%, rerunning with new random noise")

    print(ws.model_state_vector)

assert not fail, "Failed to retrieve sensor polynomial offsets grid"
