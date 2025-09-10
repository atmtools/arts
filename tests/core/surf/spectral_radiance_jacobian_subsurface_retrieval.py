import matplotlib.pyplot as plt
import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

ws.disort_settings_agendaSubsurfaceSetup()

ws.frequency_grid = np.linspace(1, 2, 11)

z = np.linspace(0, -3, 6)
tf = pyarts.arts.GeodeticField3(
    grids=[z[::-1], [0], [0]],
    grid_names=["altitude", "latitude", "longitude"],
    data=np.linspace(200, 400, len(z)).reshape(len(z), 1, 1)
)

ws.surface_fieldEarth()
ws.surface_field["t"] = tf.max()

ws.subsurface_field.bottom_depth = min(z)
ws.subsurface_field['scalar absorption'] = 5
ws.subsurface_field['scalar ssa'] = 0.99
ws.subsurface_field["t"] = tf
ws.subsurface_field["t"].alt_low = "Linear"
ws.subsurface_field["t"].alt_upp = "Linear"

NQUAD = 40

ws.spectral_radiance_transform_operatorSet(option="Tb")

ws.absorption_species = []
ws.absorption_bands = {}
ws.propagation_matrix_agendaAuto()
ws.spectral_radiance_observer_agendaSet(option="EmissionNoSensor")
ws.ray_path_observer_agendaSetGeometric()
ws.atmospheric_fieldInit(toa=0.0)


@pyarts.arts_agenda(ws=ws, fix=False)
def spectral_radiance_surface_agenda(ws):
    ws.spectral_radianceSubsurfaceDisortEmissionWithJacobian(depth_profile=z)


ws.RetrievalInit()
ws.RetrievalAddSubsurface(target="t", matrix=np.diag(1e-1 * np.ones_like(z)))

pos = [0, 0, 0]
los = [120.0, 0.0]

ws.disort_fourier_mode_dimension = 1
ws.disort_legendre_polynomial_dimension = 1
ws.disort_quadrature_dimension = 10

ws.measurement_sensorSimple(pos=pos, los=los, pol="RC")

ws.RetrievalFinalizeDiagonal()

ws.measurement_vectorFromSensor()

noise = .01
ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
epp = np.random.normal(0, noise, len(ws.frequency_grid))
ws.measurement_vector += epp
ws.subsurface_field["t"].data += 20
ws.model_state_vector_aprioriFromData()
ws.measurement_vector_fitted = []
ws.model_state_vector = []
ws.measurement_jacobian = [[]]

ws.OEM(method="gn")

plt.figure(figsize=(8, 8))

z = tf.grids[0]
plt.subplot(3, 1, 1)
plt.plot(ws.model_state_vector, z, label="fitted")
plt.plot(ws.model_state_vector_apriori, z, label="apriori")
plt.plot(ws.model_state_vector_apriori-20, z, label="true")
plt.legend()

ws.measurement_averaging_kernelCalc()

plt.subplot(3, 2, 3)
plt.plot(ws.measurement_averaging_kernel, z)
plt.plot(ws.measurement_averaging_kernel @ np.ones_like(z), z, "k")

plt.subplot(3, 2, 4)
plt.plot(ws.measurement_averaging_kernel.T, z)
plt.plot(ws.measurement_averaging_kernel.T @ np.ones_like(z), z, "k")

plt.subplot(3, 1, 3)
plt.plot(ws.frequency_grid, ws.measurement_vector, label="meas")
plt.plot(ws.frequency_grid, ws.measurement_vector_fitted, label="fitted")
plt.plot(ws.frequency_grid, ws.measurement_vector - epp, "k:", label="true")
plt.legend()

plt.tight_layout()
