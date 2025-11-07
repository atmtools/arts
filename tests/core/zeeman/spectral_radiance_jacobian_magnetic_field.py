import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

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

ws.absorption_speciesSet(species=["O2-66"])
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
ws.atmospheric_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.ray_path_observer_agendaSetGeometric()

# %% Artificial Magnetic Field

uf = pyarts.arts.FieldComponent.U
vf = pyarts.arts.FieldComponent.V
wf = pyarts.arts.FieldComponent.W
mag_u = ws.atmospheric_field["mag_u"](90e3, 0, 0)
mag_v = ws.atmospheric_field["mag_v"](90e3, 0, 0)
mag_w = ws.atmospheric_field["mag_w"](90e3, 0, 0)
ws.atmospheric_field["mag_u"] = mag_u
ws.atmospheric_field["mag_v"] = mag_v
ws.atmospheric_field["mag_w"] = mag_w
mag = {uf: mag_u, vf: mag_v, wf: mag_w}

pos = [110e3, 0, 0]
los = [160.0, 0.0]
ws.measurement_sensorSimple(pos=pos, los=los)

for fc in [uf, vf, wf]:
    ws.RetrievalInit()
    ws.RetrievalAddMagneticField(
        component=str(fc), matrix=np.diag(np.ones((1)) * 1e-10)
    )
    ws.RetrievalFinalizeDiagonal()

    fail = True

    for i in range(LIMIT):
        ws.atmospheric_field["mag_" + str(fc)] = mag[fc]
        ws.measurement_vectorFromSensor()

        ws.measurement_vector_fitted = []
        ws.model_state_vector = []
        ws.measurement_jacobian = [[]]

        ws.atmospheric_field["mag_" + str(fc)] = mag[fc] + 1e-6
        ws.model_state_vector_aprioriFromData()

        ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
        ws.measurement_vector += np.random.normal(0, noise, NF)

        ws.OEM(method="gn")

        absdiff = round(abs(mag[fc] - ws.model_state_vector[0]) * 1e9)

        print(
            f"{fc}-component: Input {round(mag[fc]*1e9)} nT, Output {round(ws.model_state_vector[0]*1e9)} nT, AbsDiff {absdiff} nT"
        )
        if absdiff >= ATOL:
            print(f"AbsDiff not less than {ATOL} nT, rerunning with new random noise")
            continue
        else:
            fail = False
            break

    assert not fail
