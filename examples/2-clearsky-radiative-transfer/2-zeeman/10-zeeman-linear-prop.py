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
        [[ 2.27853390e+02,  1.95866772e-04,  4.87491875e-05, -2.36891746e-02],
         [ 2.30930420e+02,  2.15802069e-04,  5.51979964e-05, -4.03735689e-02],
         [ 2.34870126e+02,  3.57229927e-04,  9.38280226e-05, -5.71986959e-02],
         [ 2.40419693e+02,  1.23779286e-03,  3.24557816e-04, -6.19281573e-02],
         [ 2.49835910e+02,  7.74687462e-03,  2.06857896e-03, -1.54923823e-02],
         [ 2.08948687e+02,  2.37725593e+01,  1.64969081e+00,  5.47121991e-06],
         [ 2.49835348e+02,  7.74673368e-03,  2.06869604e-03,  1.55956400e-02],
         [ 2.40418603e+02,  1.23699914e-03,  3.24403123e-04,  6.20719731e-02],
         [ 2.34868550e+02,  3.56493629e-04,  9.36665224e-05,  5.73508261e-02],
         [ 2.30928369e+02,  2.15285660e-04,  5.50826107e-05,  4.05218302e-02],
         [ 2.27850874e+02,  1.95567214e-04,  4.86829073e-05,  2.38269831e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
