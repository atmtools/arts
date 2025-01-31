import pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(118.5e9-30e9, 119e9+30e9, 14001)

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequencyByLine(fmin=118e9, fmax=120e9)
del ws.absorption_bands["O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 1 1"]

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Settings

ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.ray_path_observer_agendaSetGeometric(add_crossings=1, remove_non_crossings=1)

nlimb = 15
nup = 30
ndown = 120

# %% Checks

ws.ray_path_observersFieldProfilePseudo2D(nup=nup, ndown=ndown, nlimb=nlimb)
ws.ray_path_fieldFromObserverAgenda()
ws.spectral_flux_profileFromPathField(
    altitude_grid=ws.atmospheric_field["t"].data.grids[0]
)

ws.ray_path_atmospheric_point = [
    ws.atmospheric_field(x, 0, 0) for x in ws.atmospheric_field["t"].data.grids[0]
]

ws.nlte_line_flux_profileIntegrate()

# We are only close 
