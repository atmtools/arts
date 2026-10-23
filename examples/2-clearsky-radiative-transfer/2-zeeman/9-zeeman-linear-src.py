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
ws.rte_option = "lintau"
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
        [[ 2.27693773e+02,  4.24005564e-04,  1.10478661e-04,  5.60428826e-02],
         [ 2.30768815e+02,  6.57107391e-04,  1.71389041e-04,  6.93794910e-02],
         [ 2.34709066e+02,  1.16350287e-03,  3.03962018e-04,  9.20568215e-02],
         [ 2.40261432e+02,  2.60503017e-03,  6.82692630e-04,  1.38053297e-01],
         [ 2.49682173e+02,  9.71267320e-03,  2.56943414e-03,  2.66179941e-01],
         [ 2.09817858e+02,  2.41315155e+01,  2.23227562e+00,  5.40846395e-06],
         [ 2.49681605e+02,  9.71620266e-03,  2.57039679e-03, -2.66256405e-01],
         [ 2.40260340e+02,  2.60682925e-03,  6.83171375e-04, -1.38128865e-01],
         [ 2.34707490e+02,  1.16469914e-03,  3.04277736e-04, -9.21318124e-02],
         [ 2.30766766e+02,  6.58008544e-04,  1.71625883e-04, -6.94544638e-02],
         [ 2.27691261e+02,  4.24733194e-04,  1.10669418e-04, -5.61177737e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
