import pyarts
import numpy as np
import matplotlib.pyplot as plt

PLOT = False

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

# %% Core geometry
NQuad = 40
ws.ray_pathGeometricUplooking(latitude=0.0, longitude=0.0, max_step=1000.0)
ws.ray_path_atmospheric_pointFromPath()
ws.ray_path_frequency_gridFromPath()
ws.ray_path_propagation_matrixFromPath()

# %% Disort calculations
print("DISORT Calculations")
ws.spectral_radiance_disortClearskyDisort(NQuad=NQuad, NLeg=1)

# %% Equivalent ARTS calculations
print("ARTS Calculations")
ws.ray_pathGeometric(
    pos=[100e3, 0, 0],
    los=[180, 0],
    max_step=1000.0,
)
ws.spectral_radianceClearskyEmission()

# %% Plot results

if PLOT:
    plt.semilogy(
        ws.frequency_grid - line_f0,
        ws.spectral_radiance_disort[:, -1, (NQuad // 2) :],
        label="disort",
    )
    plt.semilogy(
        ws.frequency_grid - line_f0, ws.spectral_radiance[:, 0], "k--", lw=3
    )
    plt.semilogy(
        ws.frequency_grid - line_f0,
        ws.spectral_radiance_disort[:, -1, NQuad // 2],
        "g:",
        lw=3,
    )
    plt.semilogy(
        ws.frequency_grid - line_f0,
        ws.spectral_radiance_disort[:, -1, -1],
        "m:",
        lw=3,
    )
    plt.ylabel("Spectral radiance [W sr$^{-1}$ m$^{-2}$ Hz$^{-1}$]")
    plt.xlabel("Dirac frequency [index count]")

# %% The last test should be that we are close to the correct values

assert np.allclose(
    ws.spectral_radiance_disort[:, -1, -1] / ws.spectral_radiance[:, 0],
    1, rtol=1e-3,
), "Bad results, clearsky calculations are not close between DISORT and ARTS"
