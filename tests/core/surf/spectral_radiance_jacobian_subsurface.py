import matplotlib.pyplot as plt
import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

ws.disort_settings_agendaSubsurfaceSetup()

ws.freq_grid = [100e9]

z = np.linspace(0, -30, 401)
tf = pyarts.arts.GeodeticField3(
    grids=[z[::-1], [0], [0]],
    grid_names=["altitude", "latitude", "longitude"],
    data=np.linspace(200, 400, len(z)).reshape(len(z), 1, 1)
)

ws.surf_fieldEarth()

ws.subsurf_field.bottom_depth = min(z)
ws.subsurf_field['scalar absorption'] = 3
ws.subsurf_field['scalar ssa'] = 0.99
ws.subsurf_field["t"] = tf
ws.subsurf_field["t"].alt_low = "Linear"
ws.subsurf_field["t"].alt_upp = "Linear"

NQUAD = 40

ws.spectral_radiance_transform_operatorSet(option="Tb")

ws.abs_species = []
ws.abs_bands = {}
ws.spectral_propmat_agendaAuto()
ws.spectral_radiance_observer_agendaSet(option="EmissionNoSensor")
ws.ray_path_observer_agendaSetGeometric()
ws.atm_fieldInit(toa=0.0)

ws.RetrievalInit()
ws.RetrievalAddSubsurface(target="t", matrix=np.diag(5 * np.ones((len(z)))))
ws.measurement_sensor = []
ws.RetrievalFinalizeDiagonal()

ws.ray_path_point.pos = [0, 0, 0]
ws.ray_path_point.los = [0, 0]

ws.spectral_radianceSubsurfaceDisortEmissionWithJacobian(disort_fourier_mode_dimension=1,
                                                         disort_legendre_polynomial_dimension=1,
                                                         disort_quadrature_dimension=NQUAD,
                                                         depth_profile=z)
ws.spectral_radianceApplyUnit()

j = ws.spectral_radiance_jacobian.flatten()[::4]
plt.semilogx(j, tf.grids[0])
print(ws.spectral_radiance)
