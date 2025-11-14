import matplotlib.pyplot as plt
import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

ws.disort_settings_agendaSubsurfaceSetup()

ws.freq_grid = [100e9]

z = np.linspace(0, -30, 101)
tf = pyarts.arts.GeodeticField3(
    grids=[z[::-1], [0], [0]],
    grid_names=["altitude", "latitude", "longitude"],
    data=np.linspace(200, 400, len(z)).reshape(len(z), 1, 1)
)

ws.surf_fieldEarth()

ws.subsurf_field.bottom_depth = min(z)
ws.subsurf_field['scalar absorption'] = 0.5
ws.subsurf_field['scalar ssa'] = 0.9
ws.subsurf_field["t"] = tf
ws.subsurf_field["t"].alt_low = "Linear"
ws.subsurf_field["t"].alt_upp = "Linear"

NQUAD = 40

ws.spectral_rad_transform_operatorSet(option="Tb")

ws.abs_species = []
ws.abs_bands = {}
ws.spectral_propmat_agendaAuto()
ws.spectral_rad_observer_agendaSet(option="EmissionNoSensor")
ws.ray_path_observer_agendaSetGeometric()
ws.atm_fieldInit(toa=0.0)



ws.ray_path_point.pos = [0, 0, 0]
ws.ray_path_point.los = [0, 0]

zas = np.linspace(0, 180, 30)
data = []
for za in zas:
    ws.ray_path_point.los[0] = za
    ws.spectral_radSubsurfaceDisortEmission(disort_fourier_mode_dimension=1,
                                                 disort_legendre_polynomial_dimension=1,
                                                 disort_quadrature_dimension=NQUAD,
                                                 depth_profile=z)
    ws.spectral_rad_jac
    ws.spectral_radApplyUnit()

    data.append(ws.spectral_rad[0, 0])

plt.plot(zas, data)

assert np.all(np.sort(data) == data[::-1])
