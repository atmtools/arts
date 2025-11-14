import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt


PLOT = False  # Plot for debug
LIMIT = 50  # Rerun limit for finding a fit
NFREQ = 1001
noise = 0.1
DMOD = [0.1, 2]
ATOLS = [0.01, 1]  # The level of error tolerance
SENS = [[0, 1], [1, 2, 3]]  # The index where we are sensitive

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.freq_grid = np.linspace(-40e6, 40e6, NFREQ) + line_f0

# %% Species and line absorption

ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.abs_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

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

ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.ray_path_observer_agendaSetGeometric()

# %% Artificial VMR

vmrs = [0.2, 0.2, 0.2]
vmr_grid = pyarts.arts.GriddedField3(
    name="VMR",
    data=np.array(vmrs).reshape(3, 1, 1),
    grid_names=["Altitude", "Latitude", "Longitude"],
    grids=[[0, 50e3, 120e3], [0], [0]],
)

ws.atm_field[pyarts.arts.SpeciesEnum.O2] = vmr_grid

# %% Artificial Temperature

temps = [295, 200, 250, 200, 300]
temp_grid = pyarts.arts.GriddedField3(
    name="Temperature",
    data=np.array(temps).reshape(5, 1, 1),
    grid_names=["Altitude", "Latitude", "Longitude"],
    grids=[[0, 10e3, 30e3, 80e3, 120e3], [0], [0]],
)

ws.atm_field[pyarts.arts.AtmKey.temperature] = temp_grid
ws.atm_field[pyarts.arts.AtmKey.temperature].lat_low = "Nearest"
ws.atm_field[pyarts.arts.AtmKey.temperature].lat_upp = "Nearest"
ws.atm_field[pyarts.arts.AtmKey.temperature].lon_low = "Nearest"
ws.atm_field[pyarts.arts.AtmKey.temperature].lon_upp = "Nearest"

# %% Set up sensor

pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.measurement_sensorSimple(pos=pos, los=los)

# %% Jacobian

ws.RetrievalInit()
ws.RetrievalAddSpeciesVMR(species="O2", matrix=np.diag(np.ones((3)) * 5))
ws.RetrievalAddTemperature(matrix=np.diag(np.ones((5)) * 10))
ws.RetrievalFinalizeDiagonal()

model_statec_vector_targets = [vmrs, temps]

# %% Retrieval agenda

@pyarts.workspace.arts_agenda(ws=ws, fix=True)
def inversion_iterate_agenda(ws):
    ws.UpdateModelStates()
    ws.measurement_vectorFromSensor()
    ws.measurement_vector_fittedFromMeasurement()

# %% Conditional helper

def condition(x, xas):
    i0 = 0
    for ind in range(len(xas)):
        xa = np.array(xas[ind])
        atol = ATOLS[ind]
        sens = SENS[ind]
        n = len(xa)
        ks = [i0 + k for k in sens]
        print ("CF", x[ks], "vs", xa[sens], "for tolerance", atol)
        if not np.allclose(x[ks], xa[sens], atol=atol):
            return False
        i0 += n
    return True


# %% Core calculations

works = False

for i in range(LIMIT):
    ws.atm_field[pyarts.arts.SpeciesEnum.O2].data = vmr_grid
    ws.atm_field[pyarts.arts.AtmKey.t].data = temp_grid
    ws.measurement_vectorFromSensor()

    ws.measurement_vector_fitted = []
    ws.model_state_vector = []
    
    # Set the relevant atmospheric state to "weird", aka much different
    ws.atm_field[pyarts.arts.SpeciesEnum.O2].data += 0.1
    ws.atm_field[pyarts.arts.AtmKey.t].data += 20
    ws.model_state_vector_aprioriFromData()

    ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
    ws.measurement_vector += np.random.normal(0, noise, NFREQ)

    # Must be reset to be sure it does nothing weird inside OEM
    ws.measurement_jacobian = [[]]
    
    ws.OEM(method="gn")

    if condition(ws.model_state_vector, model_statec_vector_targets):
        works = True
        break
    print("WARNING: needed repeat run, poor condition")

if PLOT:
    f = (ws.freq_grid - line_f0) / 1e6
    plt.plot(f, ws.measurement_vector)
    plt.plot(f, ws.measurement_vector_fitted)
    plt.show()

assert works, "No more reruns allowed, failing"
