import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

# Download catalogs
pyarts.data.download()

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range
line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(-50e6, 50e6, 1001) + line_f0

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
ws.atmospheric_fieldSchmidthFieldFromIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings
ws.spectral_radiance_transform_operatorSet(option="Tb")

# %% Core calculations
pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)
ws.spectral_radianceClearskyEmission()
ws.spectral_radianceApplyUnitFromSpectralRadiance()

# %% Show results
fig, ax = plt.subplots()
ax.plot((ws.frequency_grid - line_f0) / 1e6, ws.spectral_radiance + 0)
ax.set_xlabel("Frequency offset [MHz]")
ax.set_ylabel("Spectral radiance [K]")
ax.set_title(f"Zeeman effect of {round(line_f0 / 1e6)} MHz O$_2$ line")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% Test

assert np.allclose(
    ws.spectral_radiance[::100],
    np.array(
        [
            [2.27784834e02, 4.24956142e-04, -1.07495419e-04, 5.69266632e-02],
            [2.30863853e02, 6.57653186e-04, -1.68431453e-04, 7.04138002e-02],
            [2.34806521e02, 1.16288793e-03, -3.04094356e-04, 9.33836080e-02],
            [2.40360801e02, 2.59757645e-03, -7.08915090e-04, 1.40029968e-01],
            [2.49782407e02, 9.58496203e-03, -3.00953814e-03, 2.69735122e-01],
            [2.07248447e02, 2.17886934e01, 1.95305632e00, 1.25057382e-05],
            [2.49781840e02, 9.58837415e-03, -3.01090734e-03, -2.69812575e-01],
            [2.40359708e02, 2.59935325e-03, -7.09472967e-04, -1.40106546e-01],
            [2.34804943e02, 1.16407621e-03, -3.04437307e-04, -9.34596090e-02],
            [2.30861802e02, 6.58551211e-04, -1.68679754e-04, -7.04898114e-02],
            [2.27782317e02, 4.25683062e-04, -1.07691261e-04, -5.70026253e-02],
        ]
    ),
), "Values have drifted from expected results in spectral radiance"
