import pyarts
import numpy as np
import matplotlib.pyplot as plt


PLOT = False  # Plot for debug
LIMIT = 50  # Rerun limit for finding a fit
ATOL = 0.01
NFREQ = 1001
noise = 0.1

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(-20e6, 20e6, NFREQ) + line_f0

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

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
ws.spectral_radiance_observer_agendaSet(option="Emission")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.ray_path_observer_agendaSet(option="Geometric")

# %% Artificial VMR

grid = pyarts.arts.GriddedField3(
    name="VMR",
    data=np.ones((3, 1, 1)) * 0.2,
    grid_names=["Altitude", "Latitude", "Longitude"],
    grids=[[0, 50e3, 120e3], [0], [0]],
)

ws.atmospheric_field[pyarts.arts.SpeciesEnum.O2] = grid
ws.atmospheric_field[pyarts.arts.SpeciesEnum.O2].lat_low = "Nearest"
ws.atmospheric_field[pyarts.arts.SpeciesEnum.O2].lat_upp = "Nearest"
ws.atmospheric_field[pyarts.arts.SpeciesEnum.O2].lon_low = "Nearest"
ws.atmospheric_field[pyarts.arts.SpeciesEnum.O2].lon_upp = "Nearest"


# %% Set up sensor

pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.measurement_sensorSimple(pos=pos, los=los)

# %% Jacobian

ws.RetrievalInit()
ws.RetrievalAddSpeciesVMR(species="O2", matrix=np.diag(np.ones((3)) * 5))
ws.RetrievalFinalizeDiagonal()

# %% Core calculations

for i in range(LIMIT):
    ws.atmospheric_field[pyarts.arts.SpeciesEnum.O2].data = grid
    ws.measurement_vectorFromSensor()

    ws.measurement_vector_fitted = []
    ws.model_state_vector = []
    ws.measurement_jacobian = [[]]

    ws.atmospheric_field[pyarts.arts.SpeciesEnum.O2].data += 0.1
    ws.model_state_vector_aprioriFromData()

    ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
    ws.measurement_vector += np.random.normal(0, noise, NFREQ)

    @pyarts.workspace.arts_agenda(ws=ws, fix=True)
    def inversion_iterate_agenda(ws):
        ws.UpdateModelStates()
        ws.measurement_vectorFromSensor()
        ws.measurement_vector_fittedFromMeasurement()

    ws.OEM(method="gn")

    if np.allclose(ws.model_state_vector / 0.2 - 1, 0, atol=ATOL):
        break
    print("WARNING: needed repeat run:", ws.model_state_vector)

if PLOT:
    f = (ws.frequency_grid - line_f0) / 1e6
    plt.plot(f, ws.measurement_vector)
    plt.plot(f, ws.measurement_vector_fitted)
    plt.show()

assert np.allclose(
    ws.model_state_vector / 0.2 - 1, 0, atol=ATOL
), f"This should all be {ws.model_state_vector / 0.2} within {round(ATOL*100, 1)}% of 1"
