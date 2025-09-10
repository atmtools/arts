import matplotlib.pyplot as plt
import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

ws.disort_settings_agendaSubsurfaceSetup()

ws.frequency_grid = [100e9]

ws.ray_path_point.pos = [0, 0, 0]
ws.ray_path_point.los = [0, 0]

ws.surface_fieldEarth()
ws.surface_field['t'] = 100

z = np.linspace(0, -30, 101)
tf = pyarts.arts.GeodeticField3(
    grids=[z[::-1], [0], [0]],
    grid_names=["altitude", "latitude", "longitude"],
    data=np.linspace(200, 400, len(z)).reshape(len(z), 1, 1)
)

ws.subsurface_field.bottom_depth = min(z)
ws.subsurface_field['scalar absorption'] = .5
ws.subsurface_field['scalar ssa'] = .9
ws.subsurface_field["t"] = tf
ws.subsurface_field["t"].alt_low = "Linear"
ws.subsurface_field["t"].alt_upp = "Linear"

NQUAD = 40
ws.disort_spectral_radiance_fieldSubsurfaceProfile(
    disort_fourier_mode_dimension=1,
    disort_legendre_polynomial_dimension=1,
    disort_quadrature_dimension=NQUAD,
    depth_profile=z
)

ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.disort_spectral_radiance_fieldApplyUnit()

zas = np.linspace(0, 180, 361)
data = []
for za in zas:
    ws.ray_path_point.los = [za, 0]
    ws.spectral_radianceFromDisort()
    data.append(ws.spectral_radiance[0, 0])
data = np.array(data)
plt.plot(zas, data)

assert np.all(np.sort(data) == data[::-1])
