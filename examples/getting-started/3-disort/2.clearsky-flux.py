import pyarts
import numpy as np
import matplotlib.pyplot as plt

toa = 100e3
lat = 0
lon = 0
NQuad = 40

ws = pyarts.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.frequency_grid = [line_f0]
ws.frequency_grid = np.linspace(-20e9, 2e6, 5) + line_f0

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=toa, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")
# ws.atmospheric_field[pyarts.arts.AtmKey.t] = 300.0

# %% Checks and settings

ws.spectral_radiance_unit = "Tb"
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

# %% Core Disort calculations

ws.disort_settings_agendaSet(option="Clearsky")

ws.ray_pathGeometricDownlooking(longitude=lon, latitude=lat, max_step=40_000)
ws.ray_path_atmospheric_pointFromPath()
ws.ray_path_frequency_gridFromPath()
ws.ray_path_propagation_matrixFromPath()
ws.ray_path_pointLowestFromPath()

ws.disort_spectral_flux_fieldFromAgenda(
    disort_quadrature_dimension=NQuad,
    disort_legendre_polynomial_dimension=1,
    disort_fourier_mode_dimension=1,
)

assert np.allclose(
    ws.disort_spectral_flux_field[:, :2].flatten()
    / np.array(
        [
            5.28745377e-16,
            5.21547172e-16,
            9.57947674e-18,
            2.75440511e-15,
            2.75117169e-15,
            2.67196535e-15,
            5.99210850e-16,
            5.90686327e-16,
            9.99604279e-18,
            3.03918207e-15,
            3.03530762e-15,
            2.94570409e-15,
            7.85191844e-16,
            7.73432694e-16,
            1.03766333e-17,
            3.33784022e-15,
            3.33201652e-15,
            3.21634922e-15,
            1.47946833e-15,
            1.45789557e-15,
            1.07204266e-17,
            3.65036231e-15,
            3.63385500e-15,
            3.43519461e-15,
            3.37235564e-15,
            2.86981610e-15,
            1.10269829e-17,
            3.51341265e-15,
            2.80896362e-15,
            3.94568143e-15,
        ]
    ),
    1,
), "Mismatch from historical Disort spectral fluxes"
