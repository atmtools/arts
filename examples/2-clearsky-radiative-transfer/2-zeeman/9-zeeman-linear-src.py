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
        [[ 2.27693775e+02,  4.25360933e-04,  1.02565822e-04,  5.68301097e-02],
         [ 2.30768818e+02,  6.59206212e-04,  1.59126898e-04,  7.03542474e-02],
         [ 2.34709071e+02,  1.16721699e-03,  2.82253340e-04,  9.33505189e-02],
         [ 2.40261444e+02,  2.61334904e-03,  6.34105246e-04,  1.39994110e-01],
         [ 2.49682208e+02,  9.74388124e-03,  2.38851950e-03,  2.69924367e-01],
         [ 2.09884932e+02,  2.41682141e+01,  1.73723873e+00,  5.54512081e-06],
         [ 2.49681641e+02,  9.74742248e-03,  2.38941674e-03, -2.70001907e-01],
         [ 2.40260352e+02,  2.61515400e-03,  6.34550524e-04, -1.40070740e-01],
         [ 2.34707495e+02,  1.16841714e-03,  2.82546778e-04, -9.34265633e-02],
         [ 2.30766769e+02,  6.60110278e-04,  1.59346947e-04, -7.04302732e-02],
         [ 2.27691263e+02,  4.26090913e-04,  1.02743016e-04, -5.69060525e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
