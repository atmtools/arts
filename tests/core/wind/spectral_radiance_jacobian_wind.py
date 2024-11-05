import pyarts
import numpy as np
import matplotlib.pyplot as plt

PLOT = False  # Plot for debug
LIMIT = 50  # Rerun limit for finding a fit
ATOL = 20
NF = 1001
noise = 0.5

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(-5e6, 5e6, NF) + line_f0

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequency(fmin=118e9, fmax=119e9, by_line=1)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

qid = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"
ws.absorption_bands = {qid: ws.absorption_bands[qid]}
# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_unit = "Tb"
ws.spectral_radiance_observer_agendaSet(option="EmissionUnits")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.ray_path_observer_agendaSet(option="Geometric")

# %% Artificial Wind Field

uf = pyarts.arts.FieldComponent.U
vf = pyarts.arts.FieldComponent.V
wf = pyarts.arts.FieldComponent.W
wind_u = 10.0
wind_v = 10.0
wind_w = 10.0
ws.atmospheric_field["wind_u"] = wind_u
ws.atmospheric_field["wind_v"] = wind_v
ws.atmospheric_field["wind_w"] = wind_w
wind = {uf: wind_u, vf: wind_v, wf: wind_w}

# %% Retrieval agenda


@pyarts.workspace.arts_agenda(ws=ws, fix=True)
def inversion_iterate_agenda(ws):
    ws.UpdateModelStates()
    ws.measurement_vectorFromSensor()
    ws.measurement_vector_fittedFromMeasurement()


for fc in [uf, vf, wf]:
    ws.RetrievalInit()
    ws.RetrievalAddWindField(component=str(fc), matrix=np.diag(np.ones((1)) * 100))
    ws.RetrievalFinalizeDiagonal()

    pos = [110e3, 0, 0]
    los = [140.0, 30.0]
    fail = True

    for i in range(LIMIT):
        ws.atmospheric_field["wind_u"] = wind[fc]
        ws.atmospheric_field["wind_v"] = wind[fc]
        ws.atmospheric_field["wind_w"] = wind[fc]
        ws.measurement_sensorSimple(pos=pos, los=los)
        ws.measurement_vectorFromSensor()

        ws.measurement_vector_fitted = []
        ws.model_state_vector = []
        ws.measurement_jacobian = [[]]

        ws.atmospheric_field["wind_" + str(fc)] = wind[fc] + 100
        ws.model_state_vector_aprioriFromData()

        ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
        ws.measurement_vector += np.random.normal(0, noise, NF)

        ws.OEM(method="gn")

        absdiff = round(abs(wind[fc] - ws.model_state_vector[0]))

        print(
            f"""{fc}-component:
  Apriori:  {ws.model_state_vector_apriori[0]} m/s
  Truth:    {round(wind[fc])} m/s
  Retrieved {round(ws.model_state_vector[0])} m/s
  AbsDiff   {absdiff} m/s"""
        )
        if absdiff >= ATOL:
            print(f"AbsDiff not less than {ATOL} m/s, rerunning with new random noise")
            continue
        else:
            fail = False
            break

    assert not fail
