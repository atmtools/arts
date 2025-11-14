import pyarts3 as pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.freq_grid = np.linspace(-20e6, 20e6, 101) + line_f0

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
    alt_grid=ws.atm_field["t"].data.grids[0]
)
a = np.array(ws.spectral_flux_profile)

ws.ray_path_observersFieldProfilePseudo2D(nup=2 * nup, ndown=2 * ndown, nlimb=2 * nlimb)
ws.ray_path_fieldFromObserverAgenda()
ws.spectral_flux_profileFromPathField(
    alt_grid=ws.atm_field["t"].data.grids[0]
)
b = np.array(ws.spectral_flux_profile)

assert np.allclose(a, b, rtol=1e-2)
