import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt
from copy import copy

PLOT = False  # Plot for debug
NF = 1001
noise = 1

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

ws.freq_grid = np.linspace(10e9, 400e9, NF)

# %% Species and line absorption

ws.abs_speciesSet(species=["H2O-PWR98", "O2-PWR98"])
ws.ReadCatalogData()
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)

# %% Checks and settings

ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.ray_path_observer_agendaSetGeometric()

pos = [0e3, 0, 0]
los = [20.0, 0.0]
ws.measurement_sensorSimple(pos=pos, los=los)

RAT = 0.8
field = copy(ws.atm_field["H2O"])
fieldg = copy(ws.atm_field["H2O"])

fieldg.data /= RAT
orig = np.ones(50)

ws.RetrievalInit()
ws.RetrievalAddSpeciesVMR(species="H2O", matrix=np.diag(np.ones((50))) * 1e-5)
ws.RetrievalFinalizeDiagonal()

ws.measurement_vectorFromSensor()
meas = ws.measurement_vector * 1.0
true = 1.0 * meas

ws.measurement_vector_fitted = []
ws.model_state_vector = []
ws.measurement_jacobian = [[]]

ws.atm_field["H2O"] = fieldg
ws.jacobian_targetsToggleRelativeAtmTarget(key="H2O")
ws.model_state_vector_aprioriFromData()
ws.measurement_vectorFromSensor()
apri = ws.measurement_vector * 1.0

ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
meas += np.random.normal(0, noise, NF)
ws.measurement_vector = meas

# %% OEM
ws.OEM(method="lm", lm_ga_settings=[10, 2, 2, 100, 1, 99])
ws.model_state_vectorFromData()

if PLOT:
    plt.plot(ws.freq_grid / 1e9, meas, label="orig")
    plt.plot(ws.freq_grid / 1e9, apri, label="apriori")
    plt.plot(ws.freq_grid / 1e9, ws.measurement_vector_fitted, label="fitted")
    plt.legend()
    plt.show()
    plt.plot(
        ws.model_state_vector_apriori / RAT,
        field.data.grids[0],
        ":",
        lw=3,
        label="apriori ratio",
    )
    plt.plot(ws.model_state_vector / RAT, field.data.grids[0], label="fitted ratio")
    plt.plot(ws.model_state_vector * 0 + 1, field.data.grids[0], label="true ratio")
    plt.legend()
    plt.show()

# Just a simple test that some real convergence happens (the fit is not good but better)
print(np.std(true - apri), "vs", np.std(true - ws.measurement_vector_fitted))
assert np.std(true - apri) > 2 * np.std(true - ws.measurement_vector_fitted)
