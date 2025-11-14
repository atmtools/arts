import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

PLOT = False  # Plot for debug
LIMIT = 50  # Rerun limit for finding a fit
ATOL = 0.5
NF = 11
noise = 0.5

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.freq_grid = np.linspace(-5000e6, 5000e6, NF) + line_f0

# %% Species and line absorption

ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=118e9, fmax=119e9)
ws.abs_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

bandkey = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"
ws.abs_bands = {bandkey: ws.abs_bands[bandkey]}

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atm_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.spectral_radiance_surface_agendaSet(option="SurfaceReflectance")
ws.surface_reflectance_agendaSet(option="FlatRealFresnel")
ws.ray_path_observer_agendaSetGeometric()

# %% Artificial Surface

ts = 295.0
n2 = 1.5
ws.surface_field["t"] = ts
ws.surface_field["scalar refractive index"] = n2

pos = [110e3, 0, 0]
los = [160.0, 30.0]
ws.measurement_sensorSimple(pos=pos, los=los)

ws.RetrievalInit()
ws.RetrievalAddSurface(
    target = pyarts.arts.SurfaceKey.t, matrix=np.diag(np.ones((1)) * 100)
)
ws.RetrievalFinalizeDiagonal()

fail = True

print("Retrieving temperature")
for i in range(LIMIT):
    ws.surface_field["t"] = ts
    ws.measurement_vectorFromSensor()

    ws.measurement_vector_fitted = []
    ws.model_state_vector = []
    ws.measurement_jacobian = [[]]

    ws.surface_field["t"] = ts - 30
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

print("Temperature retrieval successful\n")

# %% Artificial Surface Reset

ws.surface_field["t"] = ts
ws.surface_field["scalar refractive index"] = n2

# %% Retrieval agenda

ws.RetrievalInit()
ws.RetrievalAddSurface(
    target = "scalar refractive index", matrix=np.diag(np.ones((1)) * 1)
)
ws.RetrievalFinalizeDiagonal()

fail = True

print("Retrieving scalar refractive index")
for i in range(LIMIT):
    ws.surface_field["scalar refractive index"] = n2
    ws.measurement_vectorFromSensor()

    ws.measurement_vector_fitted = []
    ws.model_state_vector = []
    ws.measurement_jacobian = [[]]

    ws.surface_field["scalar refractive index"] = n2 + 0.3
    ws.model_state_vector_aprioriFromData()

    ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
    ws.measurement_vector += np.random.normal(0, noise, NF)

    ws.OEM(method="gn")

    absdiff = round(abs(n2 - ws.model_state_vector[0]), 2)

    print(
        f"'scalar refractive index'-component: Input {n2}, Output {round(ws.model_state_vector[0], 2)}, AbsDiff {absdiff}"
    )
    if absdiff >= ATOL:
        print(f"AbsDiff not less than {ATOL}, rerunning with new random noise")
        continue
    else:
        fail = False
        break

assert not fail

print("Flat refractive index retrieval successful\n")

# %% Artificial Surface Reset

ws.surface_field["t"] = ts
ws.surface_field["scalar refractive index"] = n2

# %% Retrieval agenda

ws.RetrievalInit()
ws.RetrievalAddSurface(
    target = "scalar refractive index", matrix=np.diag(np.ones((1)) * 1)
)
ws.RetrievalAddSurface(
    target = "t", matrix=np.diag(np.ones((1)) * 1000)
)
ws.RetrievalFinalizeDiagonal()

fail = True

print("Retrieving scalar refractive index and temperature")
ATOL = 10
for i in range(LIMIT):
    ws.surface_field["t"] = ts
    ws.surface_field["scalar refractive index"] = n2
    ws.measurement_vectorFromSensor()

    ws.measurement_vector_fitted = []
    ws.model_state_vector = []
    ws.measurement_jacobian = [[]]

    ws.surface_field["t"] = ts + 35
    ws.surface_field["scalar refractive index"] = n2 + 0.3
    ws.model_state_vector_aprioriFromData()

    ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
    ws.measurement_vector += np.random.normal(0, noise, NF)

    ws.OEM(method="gn")

    absdiff_rs = round(abs(n2 - ws.model_state_vector[0]), 2)
    absdiff_ts = round(abs(ts - ws.model_state_vector[1]), 2)

    print(
        f"'scalar refractive index'-component: Input {n2}, Output {round(ws.model_state_vector[0], 2)}, AbsDiff {absdiff_rs}"
    )
    print(
        f"t-component: Input {ts} K, Output {round(ws.model_state_vector[1], 2)} K, AbsDiff {absdiff_ts} K"
    )
    if absdiff_rs >= ATOL or absdiff_ts >= ATOL:
        print(f"AbsDiff Reflectance not less than {ATOL} %, or AbsDiff Temperature not less than {ATOL} K, rerunning with new random noise")
        continue
    else:
        fail = False
        break

assert not fail

print("Flat refractive index and temperature retrieval successful\n")
