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
            2.65884980e-15,
            2.66123187e-15,
            2.75440515e-15,
            9.57958301e-18,
            2.08163832e-17,
            5.41861020e-16,
            2.93041798e-15,
            2.93324940e-15,
            3.03918211e-15,
            9.99616715e-18,
            2.34556308e-17,
            6.14497069e-16,
            3.19240354e-15,
            3.19640875e-15,
            3.33784028e-15,
            1.03768106e-17,
            3.05931058e-17,
            8.09137668e-16,
            3.35237338e-15,
            3.36061097e-15,
            3.65036247e-15,
            1.07209030e-17,
            6.78813350e-17,
            1.56228981e-15,
            3.37235564e-15,
            2.86991074e-15,
            3.97673151e-15,
            1.52461809e-15,
            2.80896218e-15,
            3.94568143e-15,
        ]
    ),
    1,
), "Mismatch from historical Disort spectral fluxes"
