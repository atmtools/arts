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
        [[ 2.27739530e+02,  4.26740067e-04,  1.02904248e-04,  5.69651491e-02],
         [ 2.30820521e+02,  6.60713616e-04,  1.59500586e-04,  7.04540731e-02],
         [ 2.34765080e+02,  1.16940701e-03,  2.82804954e-04,  9.34185233e-02],
         [ 2.40320870e+02,  2.61832735e-03,  6.35387788e-04,  1.40054656e-01],
         [ 2.49743758e+02,  9.75522443e-03,  2.39200930e-03,  2.69770023e-01],
         [ 2.09104573e+02,  2.38722290e+01,  1.64302854e+00,  5.42896174e-06],
         [ 2.49743191e+02,  9.75876562e-03,  2.39290721e-03, -2.69847389e-01],
         [ 2.40319778e+02,  2.62013316e-03,  6.35833406e-04, -1.40131154e-01],
         [ 2.34763504e+02,  1.17060807e-03,  2.83098669e-04, -9.34944509e-02],
         [ 2.30818472e+02,  6.61619065e-04,  1.59721001e-04, -7.05300006e-02],
         [ 2.27737016e+02,  4.27471919e-04,  1.03081920e-04, -5.70409956e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
