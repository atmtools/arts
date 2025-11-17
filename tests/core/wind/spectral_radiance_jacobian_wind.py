import pyarts3 as pyarts
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
ws.freq_grid = np.linspace(-5e6, 5e6, NF) + line_f0

# %% Species and line absorption

ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=118e9, fmax=119e9)
ws.abs_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

qid = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"
ws.abs_bands = {qid: ws.abs_bands[qid]}
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

# %% Artificial Wind Field

uf = pyarts.arts.FieldComponent.U
vf = pyarts.arts.FieldComponent.V
wf = pyarts.arts.FieldComponent.W
wind_u = 10.0
wind_v = 10.0
wind_w = 10.0
ws.atm_field["wind_u"] = wind_u
ws.atm_field["wind_v"] = wind_v
ws.atm_field["wind_w"] = wind_w
wind = {uf: wind_u, vf: wind_v, wf: wind_w}


pos = [110e3, 0, 0]
los = [140.0, 30.0]
ws.measurement_sensorSimple(pos=pos, los=los)

for fc in [uf, vf, wf]:
    ws.RetrievalInit()
    ws.RetrievalAddWindField(component=str(fc), matrix=np.diag(np.ones((1)) * 100))
    ws.RetrievalFinalizeDiagonal()
    fail = True

    for i in range(LIMIT):
        ws.atm_field["wind_u"] = wind[fc]
        ws.atm_field["wind_v"] = wind[fc]
        ws.atm_field["wind_w"] = wind[fc]
        ws.measurement_vecFromSensor()

        ws.measurement_vec_fit = []
        ws.model_state_vec = []
        ws.measurement_jac = [[]]

        ws.atm_field["wind_" + str(fc)] = wind[fc] + 100
        ws.model_state_vec_aprioriFromData()

        ws.measurement_vec_error_covmatConstant(value=noise**2)
        ws.measurement_vec += np.random.normal(0, noise, NF)

        ws.OEM(method="gn")

        absdiff = round(abs(wind[fc] - ws.model_state_vec[0]))

        print(
            f"""{fc}-component:
  Apriori:  {ws.model_state_vec_apriori[0]} m/s
  Truth:    {round(wind[fc])} m/s
  Retrieved {round(ws.model_state_vec[0])} m/s
  AbsDiff   {absdiff} m/s"""
        )
        if absdiff >= ATOL:
            print(f"AbsDiff not less than {ATOL} m/s, rerunning with new random noise")
            continue
        else:
            fail = False
            break

    assert not fail
print("done!")
