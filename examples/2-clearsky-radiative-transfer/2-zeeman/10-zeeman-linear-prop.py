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
        [[ 2.27655097e+02,  4.25941643e-04,  1.02708595e-04,  5.68661397e-02],
         [ 2.30731823e+02,  6.60022963e-04,  1.59329204e-04,  7.03869074e-02],
         [ 2.34673691e+02,  1.16858889e-03,  2.82596751e-04,  9.33841959e-02],
         [ 2.40228067e+02,  2.61655084e-03,  6.34919268e-04,  1.40046082e-01],
         [ 2.49653926e+02,  9.76155919e-03,  2.39312012e-03,  2.70138698e-01],
         [ 2.09894520e+02,  2.41722890e+01,  1.73860327e+00,  5.52740752e-06],
         [ 2.49653359e+02,  9.76510257e-03,  2.39401813e-03, -2.70216174e-01],
         [ 2.40226976e+02,  2.61835630e-03,  6.35364732e-04, -1.40122633e-01],
         [ 2.34672116e+02,  1.16978940e-03,  2.82890308e-04, -9.34601655e-02],
         [ 2.30729776e+02,  6.60927425e-04,  1.59549365e-04, -7.04628681e-02],
         [ 2.27652586e+02,  4.26672085e-04,  1.02885912e-04, -5.69420287e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
