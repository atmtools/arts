import pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.Workspace()

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range
line_f0 = 118750348044.712
ws.frequency_grid = [line_f0]
ws.frequency_grid = np.linspace(-2e6, 2e6, 101) + line_f0


# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")
# ws.atmospheric_field[pyarts.arts.AtmKey.t] = 300.0

# %% Checks and settings

ws.spectral_radiance_unit = "Tb"
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

# %% Core calculations
NQuad = 40
ws.ray_pathGeometricUplooking(latitude=0.0, longitude=0.0, max_step=1000.0)
ws.ray_path_atmospheric_pointFromPath()
ws.ray_path_frequency_gridFromPath()
ws.ray_path_propagation_matrixFromPath()

ten = pyarts.arts.Tensor3()
mu = pyarts.arts.Vector()
w = pyarts.arts.Vector()
ws.spectral_radiance_disortClearskyDisort(ten, mu, w, NQuad=NQuad, NLeg=1)

ws.ray_pathGeometric(pos=[101e3, 0, 0], los=[180, 0], max_step=1000.0)
ws.spectral_radianceClearskyEmission()

plt.semilogy(
    ws.frequency_grid - line_f0, ws.spectral_radiance[:, 0], "k--", lw=5
)
plt.semilogy(
    ws.frequency_grid - line_f0, ten[:, -1, (NQuad // 2) :], label="disort"
)
plt.semilogy(ws.frequency_grid - line_f0, ten[:, -1, NQuad // 2], "g:", lw=3)
plt.semilogy(ws.frequency_grid - line_f0, ten[:, -1, -1], "m:", lw=3)
plt.ylabel("Spectral radiance [W sr$^{-1}$ m$^{-2}$ Hz$^{-1}$]")
plt.xlabel("Dirac frequency [index count]")
