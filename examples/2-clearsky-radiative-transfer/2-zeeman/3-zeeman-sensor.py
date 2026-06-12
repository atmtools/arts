import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

# Download catalogs
pyarts.data.download()

ws = pyarts.workspace.Workspace()
ws.rte_option = "constant"

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
ws.atm_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings
ws.spectral_rad_transform_operator = "Tb"
ws.ray_path_observer_agendaSetGeometric()

# %% Set up a sensor with Gaussian standard deviation channel widths on individual frequency ranges
pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.measurement_sensorSimpleGaussian(std=1e5, pos=pos, los=los, pol="RC")

# %% Core calculations
ws.measurement_vecFromSensor()

# %% Show results
fig, ax = pyarts.plot(ws.measurement_vec, xgrid=(
    ws.freq_grid - line_f0) / 1e6)
ax.set_xlabel("Frequency offset [MHz]")
ax.set_ylabel("Spectral radiance [K]")
ax.set_title(
    f"Zeeman effect of {round(line_f0 / 1e6)} MHz O$_2$ line with Gaussian channels on individual grids"
)

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% Test
assert np.allclose(
    ws.measurement_vec[::100],
    np.array(
        [227.72265043, 230.79931459, 234.76518855, 240.36673613,
         249.92055137, 211.84612808, 249.3798804 , 240.08547841,
         234.5767321 , 230.65637527, 227.60619589]
    ),
)
