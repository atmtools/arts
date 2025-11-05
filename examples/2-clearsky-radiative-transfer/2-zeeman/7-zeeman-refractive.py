import os

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

# Download catalogs
pyarts.data.download()

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 53.0669e9
ws.frequency_grid = np.linspace(-15e6, 15e6, 51) + line_f0

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=52e9, fmax=54e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/subarctic-winter/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF()

# %% Checks and settings

ws.spectral_radiance_transform_operatorSet(option="Tb")


@pyarts.arts_agenda(ws=ws, fix=True)
def propagation_matrix_single_agenda(ws):
    ws.propagation_matrix_singleInit()
    ws.propagation_matrix_singleAddVoigtLTE()


# %% Calculate and compare refractive and geometric ray paths
pos = [3571, 46, 7]
los = [20, 90]

res = []
ws.spectral_radiance_observer_position = pos
ws.spectral_radiance_observer_line_of_sight = los
ws.max_stepsize = 100.0

# %% Show results
ws.ray_path_point_back_propagation_agendaSet(option="GeometricStepwise")
ws.spectral_radianceClearskyEmissionFrequencyDependentPropagation(max_tau=1e-3)
ws.spectral_radianceApplyUnitFromSpectralRadiance(ray_path=ws.ray_paths[0])
geometric = ws.spectral_radiance * 1.0
ws.ray_path_point_back_propagation_agendaSet(option="RefractiveStepwise")
ws.spectral_radianceClearskyEmissionFrequencyDependentPropagation(max_tau=1e-3)
ws.spectral_radianceApplyUnitFromSpectralRadiance(ray_path=ws.ray_paths[0])
refractive = ws.spectral_radiance * 1.0

# %% Show results

if "ARTS_HEADLESS" not in os.environ:
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))
    freqs = ws.frequency_grid / 1e9
    pyarts.plots.StokvecVector.plot(
        geometric - refractive, fig=fig, ax=ax, freqs=freqs)
    for a in ax.flatten():
        a.legend()
        a.set_xlabel("Frequency offset [MHz]")
        a.set_ylabel("Spectral radiance [K]")
    ax[0, 0].set_title("Stokes I")
    ax[0, 1].set_title("Stokes Q")
    ax[1, 0].set_title("Stokes U")
    ax[1, 1].set_title("Stokes V")
    fig.suptitle(
        f"Difference between geometric and refractive paths around the {round(line_f0 / 1e6)} MHz O$_2$ line")
    fig.tight_layout()
    plt.show()

# %% Tests

assert np.allclose(geometric.flatten()[::21], [1.05936097e+02, -2.94133265e-04, -1.43256518e-03,  3.74373053e-01,
                                               1.14603826e+02,  9.71356380e-02, -1.35513992e-02, -3.08284667e-01,
                                               1.08782585e+02, -2.30422420e-04])

assert np.allclose(refractive.flatten()[::21], [1.05990620e+02, -2.94006494e-04, -1.43197179e-03,  3.74215042e-01,
                                                1.14655448e+02,  9.70970512e-02, -1.35451379e-02, -3.08153979e-01,
                                                1.08837328e+02, -2.30330763e-04])

assert not np.allclose(geometric, refractive)
