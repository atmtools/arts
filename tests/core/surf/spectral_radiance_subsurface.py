import matplotlib.pyplot as plt
import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

ws.disort_settings_agendaSubsurfaceSetup()

ws.frequency_grid = [100e9]

z = np.linspace(0, -30, 101)
tf = pyarts.arts.GeodeticField3(
    grids=[z[::-1], [0], [0]],
    grid_names=["altitude", "latitude", "longitude"],
    data=np.linspace(200, 400, len(z)).reshape(len(z), 1, 1)
)

ws.surface_fieldEarth()

ws.subsurface_field.bottom_depth = min(z)
ws.subsurface_field['scalar absorption'] = 0.5
ws.subsurface_field['scalar ssa'] = 0.9
ws.subsurface_field["t"] = tf
ws.subsurface_field["t"].alt_low = "Linear"
ws.subsurface_field["t"].alt_upp = "Linear"

NQUAD = 40

ws.spectral_radiance_transform_operatorSet(option="Tb")

ws.abs_species = []
ws.abs_bands = {}
ws.propagation_matrix_agendaAuto()
ws.spectral_radiance_observer_agendaSet(option="EmissionNoSensor")
ws.ray_path_observer_agendaSetGeometric()
ws.atmospheric_fieldInit(toa=0.0)



ws.ray_path_point.pos = [0, 0, 0]
ws.ray_path_point.los = [0, 0]

zas = np.linspace(0, 180, 30)
data = []
for za in zas:
    ws.ray_path_point.los[0] = za
    ws.spectral_radianceSubsurfaceDisortEmission(disort_fourier_mode_dimension=1,
                                                 disort_legendre_polynomial_dimension=1,
                                                 disort_quadrature_dimension=NQUAD,
                                                 depth_profile=z)
    ws.spectral_radiance_jacobian
    ws.spectral_radianceApplyUnit()

    data.append(ws.spectral_radiance[0, 0])

plt.plot(zas, data)

assert np.all(np.sort(data) == data[::-1])
