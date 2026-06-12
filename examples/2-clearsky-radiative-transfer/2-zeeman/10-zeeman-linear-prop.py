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
        [[ 2.27720431e+02,  1.92860863e-04,  4.80226318e-05, -2.39944534e-02],
         [ 2.30795788e+02,  2.05916096e-04,  5.28115272e-05, -4.11877682e-02],
         [ 2.34732775e+02,  6.07407437e-04,  1.53649762e-04, -3.35111787e-02],
         [ 2.40285632e+02,  1.20481776e-03,  3.16497792e-04, -6.37947118e-02],
         [ 2.49706649e+02,  5.90915798e-03,  1.64151310e-03, -9.03484330e-02],
         [ 2.09884932e+02,  2.41682141e+01,  1.73723873e+00,  5.54512226e-06],
         [ 2.49706087e+02,  5.90803897e-03,  1.64139792e-03,  9.04789486e-02],
         [ 2.40284543e+02,  1.20399960e-03,  3.16336925e-04,  6.39391756e-02],
         [ 2.34731202e+02,  6.06980643e-04,  1.53562668e-04,  3.36445494e-02],
         [ 2.30793741e+02,  2.05389194e-04,  5.26935796e-05,  4.13367589e-02],
         [ 2.27717921e+02,  1.92557235e-04,  4.79553616e-05,  2.41329191e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
