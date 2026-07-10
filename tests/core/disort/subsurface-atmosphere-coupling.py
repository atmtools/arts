import numpy as np
import pyarts3 as pa
import matplotlib.pyplot as plt

ws = pa.workspace.Workspace()

# %% Sampled frequency range
line_f0 = 118750348044.712
ws.freq_grid = np.linspace(-20e9, 2e6, 101) + line_f0

ws.abs_speciesSet(species=["H2O-PWR98", "O2-PWR98"])
ws.ReadCatalogData()

ws.spectral_propmat_agendaAuto()

ws.surf_fieldPlanet(option="Earth")
ws.surf_field['t'] = 296

z = np.linspace(0, -30, 11)
tf = pa.arts.GeodeticField3(
    grids=[z[::-1], [0], [0]],
    grid_names=["altitude", "latitude", "longitude"],
    data=np.linspace(296, 277, len(z)).reshape(len(z), 1, 1)
)
ws.subsurf_field.bottom_depth = -30
ws.subsurf_field['scalar absorption'] = 0.5
ws.subsurf_field['scalar ssa'] = 0.5
ws.subsurf_field["t"] = tf
ws.subsurf_field["t"].alt_low = "Linear"
ws.subsurf_field["t"].alt_upp = "Linear"

ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)

ws.disort_settings_agendaSubsurfaceSetup(sun_setting="None", min_optical_depth=1e-5)
ws.subsurf_disort_settings_agenda = ws.disort_settings_agenda.__copy__()
ws.disort_settings_agendaSetup()
ws.atm_disort_settings_agenda = ws.disort_settings_agenda.__copy__()

ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_rad_surface_agendaSet(option="Blackbody")

ws.disort_fourier_mode_dimension = 1
ws.disort_legendre_polynomial_dimension = 1
ws.disort_quadrature_dimension = 10
ws.alt_grid = np.append(z[::-1], np.linspace(0, 100e3, 101)[1:])

ws.disort_spectral_rad_fieldCoupledProfiles(
    tolerance=1e-8,
    max_iterations=100,
    relaxation=0.5)
ws.spectral_rad_transform_operatorSet(option="Tb")
ws.ray_point
ws.disort_spectral_rad_fieldApplyUnit()

# %% Plot results
f = ws.disort_spectral_rad_field.freq_grid / 1e9

plt.plot(f, ws.disort_spectral_rad_field.data[:, :, 0, 0])

plt.xlabel("Frequency [GHz]")
plt.ylabel("Brightness temperature [K]")
