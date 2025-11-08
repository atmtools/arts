import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

# Download catalogs
pyarts.data.download()

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range
line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(-50e6, 50e6, 1001) + line_f0

# %% Species and line absorption
ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.abs_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet
ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atm_fieldSchmidthFieldFromIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings
ws.spectral_radiance_transform_operatorSet(option="Tb")

# %% Core calculations
pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)
ws.spectral_radianceClearskyEmission()
ws.spectral_radianceApplyUnitFromSpectralRadiance()

# %% Show results
fig, ax = pyarts.plot(ws.spectral_radiance, freqs=(
    ws.frequency_grid - line_f0) / 1e6)
[a.set_xlabel("Frequency offset [MHz]") for a in ax.flatten()]
[a.set_ylabel("Spectral radiance [K]") for a in ax.flatten()]
fig.suptitle(f"Zeeman effect of {round(line_f0 / 1e6)} MHz O$_2$ line")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% Test

assert np.allclose(
    ws.spectral_radiance[::100],
    np.array(
        [[2.27784836e+02,  4.26114598e-04,  1.02751718e-04,  5.69266704e-02],
         [2.30863855e+02,  6.59892211e-04,  1.59299017e-04,  7.04138026e-02],
         [2.34806524e+02,  1.16821510e-03,  2.82508836e-04,  9.33835996e-02],
         [2.40360807e+02,  2.61610579e-03,  6.34821344e-04,  1.40029927e-01],
         [2.49782423e+02,  9.74658804e-03,  2.38967070e-03,  2.69735134e-01],
         [2.09027223e+02,  2.38234847e+01,  1.64634819e+00,  5.44981208e-06],
         [2.49781856e+02,  9.75012917e-03,  2.39056843e-03, -2.69812587e-01],
         [2.40359714e+02,  2.61791138e-03,  6.35266868e-04, -1.40106505e-01],
         [2.34804947e+02,  1.16941576e-03,  2.82802436e-04, -9.34596007e-02],
         [2.30861804e+02,  6.60797190e-04,  1.59519308e-04, -7.04898138e-02],
         [2.27782319e+02,  4.26846028e-04,  1.02929281e-04, -5.70026325e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
