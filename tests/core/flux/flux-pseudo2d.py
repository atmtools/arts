import pyarts3 as pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.frequency_grid = [line_f0]

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
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

ws.atmospheric_profileFromGrid()
ws.zenith_gridProfilePseudo2D(dza=5)
ws.spectral_radiance_fieldProfilePseudo2D()
ws.spectral_flux_profileFromSpectralRadianceField()

assert np.allclose(
    ws.spectral_flux_profile.flatten(),
    [
        7.88140316e-15,
        7.81848154e-15,
        7.70881543e-15,
        7.58215001e-15,
        7.43555846e-15,
        7.27236382e-15,
        7.10165016e-15,
        6.92661450e-15,
        6.74853314e-15,
        6.56861288e-15,
        6.38744465e-15,
        6.20570114e-15,
        6.02551781e-15,
        5.84712648e-15,
        5.67133054e-15,
        5.50406751e-15,
        5.37103772e-15,
        5.32265486e-15,
        5.36765304e-15,
        5.45655846e-15,
        5.55547727e-15,
        5.65546013e-15,
        5.74667318e-15,
        5.82297162e-15,
        5.89034270e-15,
        5.96319833e-15,
        6.10177679e-15,
        6.24863916e-15,
        6.39519578e-15,
        6.54226239e-15,
        6.68972163e-15,
        6.83690959e-15,
        6.98097689e-15,
        7.11383741e-15,
        7.20605547e-15,
        7.20884662e-15,
        7.06061226e-15,
        6.76172931e-15,
        6.34511145e-15,
        5.87636832e-15,
        5.39191993e-15,
        4.88250303e-15,
        4.27331560e-15,
        3.58806979e-15,
        3.03338367e-15,
        2.70380728e-15,
        2.52374090e-15,
        2.41833387e-15,
        2.35105965e-15,
        2.31006228e-15,
    ],
)
