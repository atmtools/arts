import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy as copy

PLOT = False  # Plot for debug
LIMIT = 50  # Rerun limit for finding a fit
ATOL = 200
NF = 1001
noise = 0.5

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(-5e6, 5e6, NF) + line_f0

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
ws.ray_path_observer_agendaSetGeometric()

# %% Artificial Magnetic Field

alt = 90e3
lat = 67.0
lon = 0.0
uf = pyarts.arts.FieldComponent.U
vf = pyarts.arts.FieldComponent.V
wf = pyarts.arts.FieldComponent.W
t3u = pyarts.arts.Tensor3([[[ws.atm_field["mag_u"](alt, lat, lon)]]])
t3v = pyarts.arts.Tensor3([[[ws.atm_field["mag_v"](alt, lat, lon)]]])
t3w = pyarts.arts.Tensor3([[[ws.atm_field["mag_w"](alt, lat, lon)]]])
g3 = ([alt], [lat], [lon])
mag_u = pyarts.arts.GeodeticField3(data=t3u, grids=g3)
mag_v = pyarts.arts.GeodeticField3(data=t3v, grids=g3)
mag_w = pyarts.arts.GeodeticField3(data=t3w, grids=g3)
ws.atm_field["mag_u"] = mag_u
ws.atm_field["mag_v"] = mag_v
ws.atm_field["mag_w"] = mag_w
ws.atm_field["mag_u"].set_extrapolation(extrapolation="Nearest")
ws.atm_field["mag_v"].set_extrapolation(extrapolation="Nearest")
ws.atm_field["mag_w"].set_extrapolation(extrapolation="Nearest")
ws.atm_fieldAbsoluteMagneticField()

pos = [110e3, 0, 0]
los = [160.0, 0.0]
ws.measurement_sensorSimple(pos=pos, los=los)
ws.measurement_vectorFromSensor()

ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
y = 1.0*ws.measurement_vector

orig = np.sqrt(mag_u.data[0, 0, 0]**2 + mag_v.data[0, 0, 0]**2 + mag_w.data[0, 0, 0]**2)

mag_u.data += 5e-6
mag_v.data += 5e-6
mag_w.data += 5e-6
ws.atm_field["mag_u"] = mag_u
ws.atm_field["mag_v"] = mag_v
ws.atm_field["mag_w"] = mag_w
ws.atm_field["mag_u"].set_extrapolation(extrapolation="Nearest")
ws.atm_field["mag_v"].set_extrapolation(extrapolation="Nearest")
ws.atm_field["mag_w"].set_extrapolation(extrapolation="Nearest")
ws.atm_fieldAbsoluteMagneticField()

modified = np.sqrt(mag_u.data[0, 0, 0]**2 +
                   mag_v.data[0, 0, 0]**2 + mag_w.data[0, 0, 0]**2)

print(
    f"Original Magnetic Field: {orig * 1e9:.2f} nT, Modified Magnetic Field: {modified * 1e9:.2f} nT"
)

ws.RetrievalInit()
ws.RetrievalAddMagneticField(component="u", matrix=np.diag(np.ones((1)) * 1e-10))
ws.RetrievalAddOverlappingMagneticField(matrix=np.diag(np.ones((1)) * 1e-10))
ws.RetrievalFinalizeDiagonal()

fail = True

for i in range(LIMIT):
    ws.measurement_vector = y + np.random.normal(0, noise, NF)
    ws.measurement_vector_fitted = []
    ws.model_state_vector = []
    ws.measurement_jacobian = [[]]

    ws.model_state_vector_aprioriFromData()

    ws.OEM(method="gn")

    absdiff = round(abs(orig - ws.model_state_vector[0]) * 1e9)

    print(
        f"Input {round(orig*1e9)} nT, Output {round(ws.model_state_vector[0]*1e9)} nT, AbsDiff {absdiff} nT"
    )
    if absdiff >= ATOL:
        print(f"AbsDiff not less than {ATOL} nT, rerunning with new random noise")
        continue
    else:
        fail = False
        break

assert not fail
