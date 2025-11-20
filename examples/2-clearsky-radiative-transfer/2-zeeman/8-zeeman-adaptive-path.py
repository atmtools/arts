import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts
from copy import deepcopy as copy

# Download catalogs
pyarts.data.download()

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range
line_f0 = 118750348044.712
ws.freq_grid = np.linspace(-50e6, 50e6, 1001) + line_f0

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
ws.atm_fieldSchmidthFieldFromIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings
ws.spectral_rad_transform_operatorSet(option="Tb")

# %% Core calculations
pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)
ws.ray_pointBackground()
ws.spectral_rad_bkgAgendasAtEndOfPath()
ws.atm_pathFromPath()
ws.freq_grid_pathFromPath()
ws.spectral_propmat_pathFromPath()
ws.spectral_radSetToBackground()
ws.spectral_radSinglePathEmissionFrequencyLoop()
ws.spectral_radApplyUnitFromSpectralRadiance()

srad0 = copy(ws.spectral_rad)
path0 = copy(ws.ray_path)

ws.spectral_propmat_pathAdaptiveHalfPath(
    max_stepsize=100., max_tau=0.05, cutoff_tau=3.0)
ws.spectral_radSetToBackground()
ws.spectral_radSinglePathEmissionFrequencyLoop()
ws.spectral_radApplyUnitFromSpectralRadiance()

freqs = (ws.freq_grid - line_f0) / 1e6
f, a = pyarts.plot(srad0, freqs=freqs, label="Regular path")
[a.set_xlabel("Frequency offset [MHz]") for a in a.flatten()]
[a.set_ylabel("Spectral radiance [K]") for a in a.flatten()]
f.suptitle(f"Zeeman effect of {round(line_f0 / 1e6)} MHz O$_2$ line")
pyarts.plot(ws.spectral_rad, freqs=freqs, fig=f, ax=a, label="Adaptive path")
[a.legend() for a in a.flatten()]

fig, ax = plt.subplots(1, 1)
d0 = path0.distances(ws.surf_field.ellipsoid)
a0 = [p.pos[0] / 1e3 for p in path0]
ax.plot(d0, a0[:-1], label="Regular path")
d1 = ws.ray_path.distances(ws.surf_field.ellipsoid)
a1 = [p.pos[0] / 1e3 for p in ws.ray_path]
ax.plot(d1, a1[:-1], label="Adaptive path")
ax.legend()
ax.set_xlabel("Distance along path [m]")
ax.set_ylabel("Altitude [km]")
fig.suptitle("Ray paths")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

assert len(ws.ray_path) > len(path0)  # More points in adaptive path
assert np.allclose([np.sum(d0)], [np.sum(d1)])
assert np.allclose(sorted(a0, reverse=True), a0)
assert np.allclose(sorted(a1, reverse=True), a1)
