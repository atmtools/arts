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
        [[ 2.27827140e+02,  4.25383147e-04,  1.02572526e-04,  5.68785875e-02],
         [ 2.30903679e+02,  6.58907241e-04,  1.59056068e-04,  7.03617869e-02],
         [ 2.34843814e+02,  1.16673604e-03,  2.82139468e-04,  9.33332773e-02],
         [ 2.40395699e+02,  2.61325506e-03,  6.34092170e-04,  1.39982571e-01],
         [ 2.49814592e+02,  9.73578709e-03,  2.38675820e-03,  2.69659250e-01],
         [ 2.08948687e+02,  2.37725593e+01,  1.64969081e+00,  5.47121870e-06],
         [ 2.49814025e+02,  9.73932755e-03,  2.38765555e-03, -2.69736783e-01],
         [ 2.40394605e+02,  2.61506011e-03,  6.34537500e-04, -1.40059220e-01],
         [ 2.34842236e+02,  1.16793609e-03,  2.82432888e-04, -9.34093435e-02],
         [ 2.30901627e+02,  6.59811580e-04,  1.59276187e-04, -7.04378734e-02],
         [ 2.27824621e+02,  4.26114017e-04,  1.02749941e-04, -5.69546563e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
