import pyarts
import numpy as np
import matplotlib.pyplot as plt

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

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_unit = "1"
ws.spectral_radiance_space_agendaSet(option="Transmission")
ws.spectral_radiance_surface_agendaSet(option="Transmission")

# %% Core calculations

pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)
ws.spectral_radianceClearskyTransmission()

# %% Show results

plt.semilogy((ws.frequency_grid - line_f0) / 1e6, ws.spectral_radiance)
plt.xlabel("Frequency offset [MHz]")
plt.ylabel("Spectral radiance [K]")
plt.title(f"Zeeman effect of {round(line_f0/1e6)} MHz O$_2$ line")

# %% Test

assert np.allclose(
    ws.spectral_radiance[::100],
    np.array(
        [
            [3.48074339e-06, 7.79626187e-11, 1.29462746e-10, -4.58316733e-08],
            [1.71316764e-06, 5.95493071e-11, 9.78633101e-11, -2.77698065e-08],
            [6.97243885e-07, 4.28547196e-11, 6.92129313e-11, -1.47794484e-08],
            [2.01811291e-07, 2.79796924e-11, 4.36322508e-11, -6.24892024e-09],
            [2.58887859e-08, 1.50722487e-11, 2.10895501e-11, -1.53914529e-09],
            [3.25340224e-13, 4.13333797e-14, -3.14840689e-13, 8.41284772e-18],
            [2.59107134e-08, 1.50971400e-11, 2.11231675e-11, 1.54394772e-09],
            [2.02108653e-07, 2.80648334e-11, 4.37633416e-11, 6.28221755e-09],
            [6.98656898e-07, 4.30420221e-11, 6.95130288e-11, 1.48859032e-08],
            [1.71750402e-06, 5.98858575e-11, 9.84130185e-11, 2.80160699e-08],
            [3.49118434e-06, 7.85004414e-11, 1.30351250e-10, 4.63074542e-08],
        ]
    ),
), "Values have drifted from expected results in spectral radiance"
