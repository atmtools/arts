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
        [[ 2.27693776e+02,  4.25360914e-04,  1.02565807e-04,  5.68301063e-02],
         [ 2.30768819e+02,  6.59206165e-04,  1.59126872e-04,  7.03542409e-02],
         [ 2.34709073e+02,  1.16721685e-03,  2.82253283e-04,  9.33505045e-02],
         [ 2.40261448e+02,  2.61334833e-03,  6.34105014e-04,  1.39994065e-01],
         [ 2.49682223e+02,  9.74387048e-03,  2.38851618e-03,  2.69924023e-01],
         [ 2.09873451e+02,  2.41579480e+01,  1.73827372e+00,  5.54548559e-06],
         [ 2.49681656e+02,  9.74741170e-03,  2.38941342e-03, -2.70001562e-01],
         [ 2.40260356e+02,  2.61515328e-03,  6.34550292e-04, -1.40070694e-01],
         [ 2.34707497e+02,  1.16841699e-03,  2.82546721e-04, -9.34265488e-02],
         [ 2.30766770e+02,  6.60110231e-04,  1.59346921e-04, -7.04302667e-02],
         [ 2.27691263e+02,  4.26090893e-04,  1.02743001e-04, -5.69060490e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
