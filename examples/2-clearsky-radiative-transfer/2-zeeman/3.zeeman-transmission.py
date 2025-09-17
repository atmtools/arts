import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

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
ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)
ws.spectral_radianceClearskyTransmission()

# %% Show results
fig, ax = plt.subplots()
ax.semilogy((ws.frequency_grid - line_f0) / 1e6, ws.spectral_radiance)
ax.set_xlabel("Frequency offset [MHz]")
ax.set_ylabel("Spectral radiance [K]")
ax.set_title(f"Zeeman effect of {round(line_f0 / 1e6)} MHz O$_2$ line")

# %% Test
assert np.allclose(
    ws.spectral_radiance[::100],
    np.array(
        [
            [3.48823498e-06, -1.46551713e-10, 3.81072431e-11, -4.59221798e-08],
            [1.71706875e-06, -1.10970204e-10, 2.94013625e-11, -2.78281718e-08],
            [6.98940567e-07, -7.87051453e-11, 2.15072772e-11, -1.48128572e-08],
            [2.02345215e-07, -4.99009177e-11, 1.44875960e-11, -6.26439080e-09],
            [2.59661967e-08, -2.45683770e-11, 8.49367757e-12, -1.54349235e-09],
            [3.26929916e-13, 2.91596383e-13, 1.29579122e-13, 8.52093039e-18],
            [2.59895225e-08, -2.46091531e-11, 8.50852019e-12, 1.54840970e-09],
            [2.02661063e-07, -5.00559190e-11, 1.45335398e-11, 6.29845922e-09],
            [7.00439370e-07, -7.90569785e-11, 2.16047322e-11, 1.49217232e-08],
            [1.72166340e-06, -1.11611698e-10, 2.95731412e-11, 2.80799246e-08],
            [3.49928820e-06, -1.47585397e-10, 3.83784757e-11, 4.64084437e-08],
        ]
    ),
), "Values have drifted from expected results in spectral radiance"
