import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

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
ws.rte_option = "constant"
ws.spectral_radClearskyEmission()
ws.spectral_radApplyUnitFromSpectralRadiance()

# %% Show results
fig, ax = pyarts.plot(ws.spectral_rad, freqs=(
    ws.freq_grid - line_f0) / 1e6)
[a.set_xlabel("Frequency offset [MHz]") for a in ax.flatten()]
[a.set_ylabel("Spectral radiance [K]") for a in ax.flatten()]
fig.suptitle(f"Zeeman effect of {round(line_f0 / 1e6)} MHz O$_2$ line")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% Test

assert np.allclose(
    ws.spectral_rad[::100],
    np.array(
        [[ 2.27651280e+02,  4.24744310e-04,  1.10673951e-04,  5.60913683e-02],
         [ 2.30728862e+02,  6.58098870e-04,  1.71652787e-04,  6.94319831e-02],
         [ 2.34671705e+02,  1.16498417e-03,  3.04360612e-04,  9.21073412e-02],
         [ 2.40226508e+02,  2.60786693e-03,  6.83473256e-04,  1.38100104e-01],
         [ 2.49649931e+02,  9.72331946e-03,  2.57251799e-03,  2.66252039e-01],
         [ 2.09846221e+02,  2.41330397e+01,  2.23159437e+00,  5.38804694e-06],
         [ 2.49649365e+02,  9.72684954e-03,  2.57348102e-03, -2.66328423e-01],
         [ 2.40225417e+02,  2.60966656e-03,  6.83952202e-04, -1.38175601e-01],
         [ 2.34670130e+02,  1.16618107e-03,  3.04676523e-04, -9.21822679e-02],
         [ 2.30726815e+02,  6.59000674e-04,  1.71889816e-04, -6.95068813e-02],
         [ 2.27648769e+02,  4.25472509e-04,  1.10864870e-04, -5.61661527e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
