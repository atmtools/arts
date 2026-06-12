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
        [[ 2.27651281e+02,  4.26102083e-04,  1.02747359e-04,  5.68792774e-02],
         [ 2.30728865e+02,  6.60200934e-04,  1.59372202e-04,  7.04074781e-02],
         [ 2.34671710e+02,  1.16870319e-03,  2.82624430e-04,  9.34017504e-02],
         [ 2.40226519e+02,  2.61619541e-03,  6.34833393e-04,  1.40041579e-01],
         [ 2.49649966e+02,  9.75456580e-03,  2.39140821e-03,  2.69997496e-01],
         [ 2.09913635e+02,  2.41699652e+01,  1.73644537e+00,  5.52423749e-06],
         [ 2.49649399e+02,  9.75810766e-03,  2.39230583e-03, -2.70074955e-01],
         [ 2.40225428e+02,  2.61800091e-03,  6.35278864e-04, -1.40118138e-01],
         [ 2.34670135e+02,  1.16990396e-03,  2.82918050e-04, -9.34777298e-02],
         [ 2.30726817e+02,  6.61105654e-04,  1.59592427e-04, -7.04834283e-02],
         [ 2.27648771e+02,  4.26832634e-04,  1.02924704e-04, -5.69551120e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
