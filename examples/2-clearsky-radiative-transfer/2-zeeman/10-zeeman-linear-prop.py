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
ws.spectral_radClearskyLinearInTauAndPropEmission()
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
        [[ 2.27830141e+02,  4.25489128e-04,  1.02599187e-04,  5.68881918e-02],
         [ 2.30907190e+02,  6.59070807e-04,  1.59097448e-04,  7.03735321e-02],
         [ 2.34847968e+02,  1.16702319e-03,  2.82212718e-04,  9.33486760e-02],
         [ 2.40400744e+02,  2.61388422e-03,  6.34254900e-04,  1.40005198e-01],
         [ 2.49821088e+02,  9.73795166e-03,  2.38733210e-03,  2.69700244e-01],
         [ 2.08949873e+02,  2.37747404e+01,  1.64966784e+00,  5.47066239e-06],
         [ 2.49820521e+02,  9.74149272e-03,  2.38822964e-03, -2.69777785e-01],
         [ 2.40399650e+02,  2.61568960e-03,  6.34700330e-04, -1.40081856e-01],
         [ 2.34846390e+02,  1.16822346e-03,  2.82506203e-04, -9.34247505e-02],
         [ 2.30905136e+02,  6.59975316e-04,  1.59317615e-04, -7.04496271e-02],
         [ 2.27827622e+02,  4.26220137e-04,  1.02776641e-04, -5.69642694e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
