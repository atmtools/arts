import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

PLOT = False  # Plot for debug
LIMIT = 50  # Rerun limit for finding a fit
ATOL = 5
NF = 11
noise = 0.5

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.freq_grid = np.linspace(-500e6, 500e6, NF) + line_f0

# %% Species and line absorption

ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=118e9, fmax=119e9)
ws.abs_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

bandkey = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"
ws.abs_bands = {bandkey: ws.abs_bands[bandkey]}

# %% Use the automatic agenda setter for propagation matrix calculations
ws.spectral_propmat_agendaAuto()

# %% Grids and planet

ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atm_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_rad_transform_operatorSet(option="Tb")
ws.ray_path_observer_agendaSetGeometric()

# %% Artificial Surface

ts = 295.0
ws.surf_field["t"] = ts

# %% Retrieval agenda

@pyarts.workspace.arts_agenda(ws=ws, fix=True)
def inversion_iterate_agenda(ws):
    ws.UpdateModelStates()
    ws.measurement_vectorFromSensor()
    ws.measurement_vector_fittedFromMeasurement()

pos = [110e3, 0, 0]
los = [160.0, 0.0]
ws.measurement_sensorSimple(pos=pos, los=los)

ws.RetrievalInit()
ws.RetrievalAddSurface(
    target = pyarts.arts.SurfaceKey.t, matrix=np.diag(np.ones((1)) * 1000)
)
ws.RetrievalFinalizeDiagonal()

fail = True

for i in range(LIMIT):
    ws.surf_field["t"] = ts
    ws.measurement_vectorFromSensor()

    ws.measurement_vector_fitted = []
    ws.model_state_vector = []
    ws.measurement_jacobian = [[]]

    ws.surf_field["t"] = ts + 30
    ws.model_state_vector_aprioriFromData()

    ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
    ws.measurement_vector += np.random.normal(0, noise, NF)

    ws.OEM(method="gn")

    absdiff = round(abs(ts - ws.model_state_vector[0]), 2)

    print(
        f"t-component: Input {ts} K, Output {round(ws.model_state_vector[0], 2)} K, AbsDiff {absdiff} K"
    )
    if absdiff >= ATOL:
        print(f"AbsDiff not less than {ATOL} K, rerunning with new random noise")
        continue
    else:
        fail = False
        break

assert not fail
