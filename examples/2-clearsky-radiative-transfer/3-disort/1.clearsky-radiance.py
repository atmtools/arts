import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

toa = 120e3
lat = 0
lon = 0
NQuad = 40

ws = pyarts.Workspace()

# %% Sampled frequency range
line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(-20e9, 2e6, 101) + line_f0

# %% Species and line absorption
ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet
ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=toa, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings
ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

# %% Core Disort calculations
ws.disort_settings_agendaSetup()

ws.disort_spectral_radiance_fieldProfile(
    longitude=lon,
    latitude=lat,
    disort_quadrature_dimension=NQuad,
    disort_legendre_polynomial_dimension=1,
    disort_fourier_mode_dimension=1,
)

# %% Equivalent ARTS calculations
ws.ray_pathGeometricDownlooking(
    latitude=lat,
    longitude=lon,
    max_step=1000.0,
)
ws.spectral_radianceClearskyEmission()

# %% Plot results
f = ws.frequency_grid - line_f0

fig, ax = plt.subplots()
ax.semilogy(
    f,
    ws.disort_spectral_radiance_field.data[:, 0, 0, : (NQuad // 2)],
    label="disort",
)
ax.semilogy(f, ws.spectral_radiance[:, 0], "k--", lw=3)
ax.semilogy(
    f,
    ws.disort_spectral_radiance_field.data[:, 0, 0, (NQuad // 2) - 1],
    "g:",
    lw=3,
)
ax.semilogy(
    f,
    ws.disort_spectral_radiance_field.data[:, 0, 0, 0],
    "m:",
    lw=3,
)
ax.set_ylabel("Spectral radiance [W sr$^{-1}$ m$^{-2}$ Hz$^{-1}$]")
ax.set_xlabel("Dirac frequency [index count]")
ax.set_title("Downlooking")

# %% The last test should be that we are close to the correct values
assert np.allclose(
    ws.disort_spectral_radiance_field.data[:, 0, 0, 0] / ws.spectral_radiance[:, 0],
    1,
    rtol=1e-3,
), "Bad results, clearsky calculations are not close between DISORT and ARTS"
