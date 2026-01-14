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
        [[ 2.27827137e+02,  4.25383355e-04,  1.02572576e-04, 5.68786105e-02],
         [ 2.30903675e+02,  6.58907710e-04,  1.59056181e-04, 7.03618286e-02],
         [ 2.34843808e+02,  1.16673733e-03,  2.82139773e-04, 9.33333640e-02],
         [ 2.40395689e+02,  2.61325977e-03,  6.34093204e-04, 1.39982788e-01],
         [ 2.49814591e+02,  9.73578780e-03,  2.38675267e-03, 2.69659401e-01],
         [ 2.08998520e+02,  2.38211765e+01,  1.64703709e+00, 5.47032781e-06],
         [ 2.49814023e+02,  9.73932829e-03,  2.38765002e-03, -2.69736936e-01],
         [ 2.40394595e+02,  2.61506483e-03,  6.34538536e-04,-1.40059437e-01],
         [ 2.34842230e+02,  1.16793738e-03,  2.82433194e-04,-9.34094305e-02],
         [ 2.30901623e+02,  6.59812051e-04,  1.59276300e-04, -7.04379152e-02],
         [ 2.27824619e+02,  4.26114226e-04,  1.02749992e-04,-5.69546795e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
