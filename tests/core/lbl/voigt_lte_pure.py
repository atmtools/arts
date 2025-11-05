import pyarts3 as pyarts
import numpy as np
import os

ws = pyarts.Workspace()

# %% Sampled frequency range
line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(-5e6, 5e6, 1001) + line_f0

# %% Species and line absorption
ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequencyByLine(fmin=40e9, fmax=1200e9)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Set the atmosphere
ws.atmospheric_point["t"] = 300.0
ws.atmospheric_point["p"] = 1013.25e-2
ws.atmospheric_point["O2"] = 0.2095
ws.atmospheric_point.mag = [10e-6, 50e-6, 1e-6]
ws.ray_path_point.los = [180.0, 0.0]

# %% Set the jacobian
ws.jacobian_targetsInit()
ws.jacobian_targetsAddTemperature()
ws.jacobian_targetsAddSpeciesVMR(species="O2")
ws.jacobian_targetsAddMagneticField(component="u")
ws.jacobian_targetsAddMagneticField(component="v")
ws.jacobian_targetsAddMagneticField(component="w")
ws.jacobian_targetsAddWindField(component="u")
ws.jacobian_targetsAddWindField(component="v")
ws.jacobian_targetsAddWindField(component="w")
x = ["Temperature", "VMR_O2", "Mag_u", "Mag_v", "Mag_w", "Wind_u", "Wind_v", "Wind_w"]

# %% Calculations

ws.propagation_matrixInit()
ws.propagation_matrixAddVoigtLTE()

d = ws.dispersion * 1.0
pm = ws.propagation_matrix * 1.0
dpm = ws.propagation_matrix_jacobian * 1.0
dd = ws.dispersion_jacobian * 1.0
f = ws.frequency_grid * 1.0
ws.frequency_grid = f + 10
ws.propagation_matrixInit()
ws.propagation_matrixAddVoigtLTE()
d2 = ws.dispersion * 1.0
pm2 = ws.propagation_matrix * 1.0
dpm2 = (pm2 - pm) / 10
ws.frequency_grid = f

# %% Standard code

ws.propagation_matrix = []
ws.propagation_matrix_jacobian = [[]]
ws.propagation_matrixInit()
ws.propagation_matrixAddLines()

assert np.allclose(pm, ws.propagation_matrix)
# assert np.allclose(dpm, ws.propagation_matrix_jacobian)  # Disabled due zero-crossings

if "ARTS_HEADLESS" not in os.environ:
  import matplotlib.pyplot as plt
  fig, ax = plt.subplots(3, 1, figsize=(8, 12))
  ax[0].plot(ws.frequency_grid/1e9, d)
  ax[1].plot(ws.frequency_grid/1e9, pm)
  ax[2].plot(ws.frequency_grid/1e9, ws.propagation_matrix )

  for i in range(dd.shape[0]):
    fig, ax = plt.subplots(3 + ("Wind" in x[i]), 1, figsize=(8, 12+ ("Wind" in x[i])*4))
    ax[0].plot(ws.frequency_grid/1e9, dd[i])
    ax[1].plot(ws.frequency_grid/1e9, dpm[i])
    if "Wind" in x[i]:
      ax[3].plot(ws.frequency_grid/1e9, dpm2)
    ax[2].plot(ws.frequency_grid/1e9, ws.propagation_matrix_jacobian[i])
    ax[0].set_title(f"Jacobian for {x[i]}")

