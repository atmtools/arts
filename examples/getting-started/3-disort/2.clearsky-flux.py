import pyarts
import numpy as np
import matplotlib.pyplot as plt

PLOT = False
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

ws.ray_pathGeometricUplooking(longitude=lon, latitude=lat, max_step=40_000)
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
    / [
        9.59633877e-18,
        9.57950338e-18,
        9.57947235e-18,
        2.68124574e-15,
        2.68124078e-15,
        2.68124078e-15,
        1.00156763e-17,
        9.99604747e-18,
        9.99601259e-18,
        2.95621931e-15,
        2.95621356e-15,
        2.95621353e-15,
        1.04047833e-17,
        1.03765939e-17,
        1.03766146e-17,
        3.22890938e-15,
        3.22890135e-15,
        3.22890140e-15,
        1.07961162e-17,
        1.07204281e-17,
        1.07204267e-17,
        3.44348175e-15,
        3.44346269e-15,
        3.44346269e-15,
        3.32971703e-15,
        1.90723245e-17,
        1.10269829e-17,
        3.45080504e-15,
        2.56537096e-15,
        2.56528602e-15,
    ],
    1,
), "Mismatch from historical Disort spectral fluxes"
