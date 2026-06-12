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
ws.freq_grid = np.linspace(-50e6, 50e6, nf) + line_f0

# %% Species and line absorption
ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.abs_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.spectral_propmat_agendaAuto()

# %% Grids and planet
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atm_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings
ws.spectral_rad_space_agendaSet(option="Transmission")
ws.spectral_rad_surface_agendaSet(option="Transmission")

# %% Core calculations
pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)
ws.spectral_radClearskyTransmission()

# %% Show results
fig, ax = pyarts.plot(ws.spectral_rad, freqs=(
    ws.freq_grid - line_f0) / 1e6, component='I')
ax.set_yscale('log')
ax.set_xlabel("Frequency offset [MHz]")
ax.set_ylabel("Spectral radiance [K]")
ax.set_title(f"Zeeman effect of {round(line_f0 / 1e6)} MHz O$_2$ line")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% Test
assert np.allclose(
    ws.spectral_rad[::100],
    np.array(
        [[ 3.55416974e-06, -1.49626896e-10, -3.59945621e-11, -4.65675860e-08],
         [ 1.75504614e-06, -1.14112664e-10, -2.74625108e-11, -2.83292389e-08],
         [ 7.16577429e-07, -8.14996162e-11, -1.96228377e-11, -1.51461032e-08],
         [ 2.07916647e-07, -5.19869800e-11, -1.25199070e-11, -6.43052931e-09],
         [ 2.67056326e-08, -2.58535066e-11, -6.20359298e-12, -1.58692924e-09],
         [ 7.87383592e-13,  7.80210615e-13, -6.46005354e-14,  8.84600141e-18],
         [ 2.67295620e-08, -2.58964360e-11, -6.21391407e-12,  1.59197636e-09],
         [ 2.08240338e-07, -5.21480127e-11, -1.25589127e-11,  6.46545336e-09],
         [ 7.18110068e-07, -8.18625471e-11, -1.97105925e-11,  1.52574035e-08],
         [ 1.75973088e-06, -1.14770061e-10, -2.76211969e-11,  2.85857722e-08],
         [ 3.56540626e-06, -1.50680035e-10, -3.62484502e-11,  4.70615011e-08]]
    ),
), "Values have drifted from expected results in spectral radiance"
