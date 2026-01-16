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
ws.rte_option = "linprop"
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
        [[2.27853388e+02,  1.95866945e-04,  4.87492271e-05, -2.36891567e-02],
         [2.30930417e+02,  2.15802433e-04,  5.51980780e-05, -4.03735389e-02],
         [2.34870122e+02,  3.57230883e-04,  9.38282269e-05, -5.71986370e-02],
         [2.40419686e+02,  1.23779635e-03,  3.24558457e-04, -6.19280114e-02],
         [2.49835914e+02,  7.74686776e-03,  2.06857040e-03, -1.54924878e-02],
         [2.08998520e+02,  2.38211765e+01,  1.64703709e+00,  5.47032904e-06],
         [2.49835352e+02,  7.74672684e-03,  2.06868749e-03,  1.55957453e-02],
         [2.40418596e+02,  1.23700264e-03,  3.24403765e-04,  6.20718269e-02],
         [2.34868546e+02,  3.56494587e-04,  9.36667274e-05,  5.73507670e-02],
         [2.30928366e+02,  2.15286025e-04,  5.50826926e-05,  4.05218001e-02],
         [2.27850872e+02,  1.95567388e-04,  4.86829472e-05,  2.38269651e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
