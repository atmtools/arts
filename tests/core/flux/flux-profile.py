import pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(-20e6, 20e6, 101) + line_f0

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

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
a = np.array(ws.spectral_flux_profile)

ws.ray_path_observersFieldProfilePseudo2D(nup=2 * nup, ndown=2 * ndown, nlimb=2 * nlimb)
ws.ray_path_fieldFromObserverAgenda()
ws.spectral_flux_profileFromPathField(
    altitude_grid=ws.atmospheric_field["t"].data.grids[0]
)
b = np.array(ws.spectral_flux_profile)

assert np.allclose(a, b, rtol=1e-2)
