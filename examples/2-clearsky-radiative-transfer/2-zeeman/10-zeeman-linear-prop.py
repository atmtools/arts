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
        [[ 2.27720432e+02,  1.92860851e-04,  4.80226218e-05, -2.39944522e-02],
         [ 2.30795789e+02,  2.05916067e-04,  5.28115105e-05, -4.11877660e-02],
         [ 2.34732777e+02,  6.07407341e-04,  1.53649729e-04, -3.35111742e-02],
         [ 2.40285634e+02,  1.20481729e-03,  3.16497679e-04, -6.37946966e-02],
         [ 2.49706659e+02,  5.90915050e-03,  1.64151117e-03, -9.03483435e-02],
         [ 2.09873451e+02,  2.41579480e+01,  1.73827372e+00,  5.54548706e-06],
         [ 2.49706097e+02,  5.90803149e-03,  1.64139599e-03,  9.04788589e-02],
         [ 2.40284545e+02,  1.20399913e-03,  3.16336811e-04,  6.39391603e-02],
         [ 2.34731203e+02,  6.06980547e-04,  1.53562635e-04,  3.36445449e-02],
         [ 2.30793742e+02,  2.05389164e-04,  5.26935628e-05,  4.13367566e-02],
         [ 2.27717922e+02,  1.92557223e-04,  4.79553515e-05,  2.41329179e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
