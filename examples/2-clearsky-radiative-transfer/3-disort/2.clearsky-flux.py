import pyarts3 as pyarts
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

ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=toa, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")
# ws.atmospheric_field[pyarts.arts.AtmKey.t] = 300.0

# %% Checks and settings

ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

# %% Core Disort calculations

ws.disort_settings_agendaSetup()

ws.ray_pathGeometricDownlooking(longitude=lon, latitude=lat, max_step=40_000)

ws.disort_spectral_flux_fieldFromAgenda(
    disort_quadrature_dimension=NQuad,
    disort_legendre_polynomial_dimension=1,
    disort_fourier_mode_dimension=1,
)

assert np.allclose(
    np.append(ws.disort_spectral_flux_field.up,
              ws.disort_spectral_flux_field.down_diffuse, axis=1).flatten()
    / np.array(
        [
            2.65924430e-15,
            2.66162733e-15,
            2.75440515e-15,
            9.57958182e-18,
            2.08076327e-17,
            5.39749082e-16,
            2.93074072e-15,
            2.93357370e-15,
            3.03918211e-15,
            9.99616600e-18,
            2.34511632e-17,
            6.12773186e-16,
            3.19264874e-15,
            3.19665668e-15,
            3.33784028e-15,
            1.03768102e-17,
            3.05943372e-17,
            8.07853130e-16,
            3.35251230e-15,
            3.36075470e-15,
            3.65036247e-15,
            1.07209025e-17,
            6.78926741e-17,
            1.56163327e-15,
            3.37235023e-15,
            2.86992821e-15,
            3.97673151e-15,
            1.52446216e-15,
            2.80896759e-15,
            3.94566396e-15,
        ]
    ),
    1,
), "Mismatch from historical Disort spectral fluxes"
