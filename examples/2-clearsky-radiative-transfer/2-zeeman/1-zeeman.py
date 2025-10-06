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
ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet
ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldSchmidthFieldFromIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings
ws.spectral_radiance_transform_operatorSet(option="Tb")

# %% Core calculations
pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)
ws.spectral_radianceClearskyEmission()
ws.spectral_radianceApplyUnitFromSpectralRadiance()

# %% Show results
fig, ax = plt.subplots()
ax.plot((ws.frequency_grid - line_f0) / 1e6, ws.spectral_radiance + 0)
ax.set_xlabel("Frequency offset [MHz]")
ax.set_ylabel("Spectral radiance [K]")
ax.set_title(f"Zeeman effect of {round(line_f0 / 1e6)} MHz O$_2$ line")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% Test

assert np.allclose(
    ws.spectral_radiance[::100],
    np.array(
        [[2.27784836e+02, -4.26290412e-04,  1.02016478e-04, -5.69266704e-02],
         [2.30863855e+02, -6.60229454e-04,  1.57888845e-04, -7.04138026e-02],
         [2.34806524e+02, -1.16900653e-03,  2.79199852e-04, -9.33835997e-02],
         [2.40360807e+02, -2.61876297e-03,  6.23712613e-04, -1.40029927e-01],
         [2.49782423e+02, -9.76705620e-03,  2.30408638e-03, -2.69735137e-01],
         [2.09029052e+02, -2.19472118e+01,  9.41746423e+00, -5.45169624e-06],
         [2.49781856e+02, -9.77061823e-03,  2.30489683e-03,  2.69812591e-01],
         [2.40359714e+02, -2.62057379e-03,  6.24136297e-04,  1.40106505e-01],
         [2.34804947e+02, -1.17020952e-03,  2.79483770e-04,  9.34596007e-02],
         [2.30861804e+02, -6.61135750e-04,  1.58103635e-04,  7.04898138e-02],
         [2.27782319e+02, -4.27022704e-04,  1.02190442e-04,  5.70026325e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
