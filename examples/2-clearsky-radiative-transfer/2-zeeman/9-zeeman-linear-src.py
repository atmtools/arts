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
ws.spectral_radClearskyLinearInTauEmission()
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
        [[2.27742531e+02,  4.26846049e-04,  1.02930910e-04, 5.69747534e-02],
         [2.30824032e+02,  6.60877182e-04,  1.59541966e-04, 7.04658183e-02],
         [2.34769234e+02,  1.16969415e-03,  2.82878204e-04, 9.34339220e-02],
         [2.40325915e+02,  2.61895652e-03,  6.35550519e-04, 1.40077284e-01],
         [2.49750254e+02,  9.75738899e-03,  2.39258319e-03, 2.69811017e-01],
         [2.09105759e+02,  2.38744101e+01,  1.64300556e+00, 5.42840543e-06],
         [2.49749687e+02,  9.76093079e-03,  2.39348131e-03, -2.69888391e-01],
         [2.40324823e+02,  2.62076265e-03,  6.35996235e-04, -1.40153790e-01],
         [2.34767658e+02,  1.17089544e-03,  2.83171984e-04, -9.35098578e-02],
         [2.30821982e+02,  6.61782801e-04,  1.59762430e-04, -7.05417543e-02],
         [2.27740017e+02,  4.27578039e-04,  1.03108621e-04, -5.70506086e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
