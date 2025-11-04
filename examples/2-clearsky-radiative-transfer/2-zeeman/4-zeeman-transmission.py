import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

# Download catalogs
pyarts.data.download()

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range
line_f0 = 118750348044.712
nf = 1001
ws.frequency_grid = np.linspace(-50e6, 50e6, nf) + line_f0

# %% Species and line absorption
ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet
ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings
ws.spectral_radiance_space_agendaSet(option="Transmission")
ws.spectral_radiance_surface_agendaSet(option="Transmission")

# %% Core calculations
pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)
ws.spectral_radianceClearskyTransmission()

# %% Show results
fig, ax = plt.subplots()
ax.semilogy((ws.frequency_grid - line_f0) / 1e6, ws.spectral_radiance)
ax.set_xlabel("Frequency offset [MHz]")
ax.set_ylabel("Spectral radiance [K]")
ax.set_title(f"Zeeman effect of {round(line_f0 / 1e6)} MHz O$_2$ line")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% Test
assert np.allclose(
    ws.spectral_radiance[::100],
    np.array(
        [[3.48824591e-06, -1.47203439e-10, -3.54155511e-11, -4.59223182e-08],
         [1.71707407e-06, -1.11587582e-10, -2.68557096e-11, -2.78282572e-08],
         [6.98942589e-07, -7.92932217e-11, -1.90895244e-11, -1.48129013e-08],
         [2.02345668e-07, -5.04729873e-11, -1.21507664e-11, -6.26440637e-09],
         [2.59661793e-08, -2.51846307e-11, -6.03516645e-12, -1.54349340e-09],
         [6.89347083e-13,  6.83290031e-13, -5.69298697e-14,  7.52268204e-18],
         [2.59895049e-08, -2.52264844e-11, -6.04518733e-12,  1.54841075e-09],
         [2.02661517e-07, -5.06298315e-11, -1.21887251e-11,  6.29847485e-09],
         [7.00441394e-07, -7.96478027e-11, -1.91752558e-11,  1.49217675e-08],
         [1.72166873e-06, -1.12232853e-10, -2.70115093e-11,  2.80800108e-08],
         [3.49929916e-06, -1.48242055e-10, -3.56660238e-11,  4.64085834e-08]]
    ),
), "Values have drifted from expected results in spectral radiance"
