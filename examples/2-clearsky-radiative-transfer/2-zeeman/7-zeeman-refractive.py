import os

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

# Download catalogs
pyarts.data.download()

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 53.0669e9
ws.freq_grid = np.linspace(-15e6, 15e6, 51) + line_f0

# %% Species and line absorption

ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.abs_bandsSetZeeman(species="O2-66", fmin=52e9, fmax=54e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.spectral_propmat_agendaAuto()

# %% Grids and planet

ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atm_fieldIGRF()

# %% Checks and settings

ws.spectral_rad_transform_operatorSet(option="Tb")


@pyarts.arts_agenda(ws=ws, fix=True)
def single_propmat_agenda(ws):
    ws.single_propmatInit()
    ws.single_propmatAddVoigtLTE()


# %% Calculate and compare refractive and geometric ray paths
pos = [3571, 46, 7]
los = [20, 90]

res = []
ws.obs_pos = pos
ws.obs_los = los
ws.max_stepsize = 20000.0

# %% Show results
ws.ray_point_back_propagation_agendaSet(option="GeometricStepwise")
ws.spectral_radClearskyEmissionFrequencyDependentPropagation(max_tau=1e-2)
ws.spectral_radApplyUnitFromSpectralRadiance(ray_path=ws.spectral_ray_path[0])
geometric = ws.spectral_rad * 1.0
ws.ray_point_back_propagation_agendaSet(option="RefractiveStepwise")
ws.spectral_radClearskyEmissionFrequencyDependentPropagation(max_tau=1e-2)
ws.spectral_radApplyUnitFromSpectralRadiance(ray_path=ws.spectral_ray_path[0])
refractive = ws.spectral_rad * 1.0

# %% Show results

if "ARTS_HEADLESS" not in os.environ:
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))
    freqs = ws.freq_grid / 1e9
    pyarts.plots.StokvecVector.plot(
        geometric - refractive, fig=fig, ax=ax, freqs=freqs)
    for a in ax.flatten():
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

assert np.allclose(geometric.flatten()[::21], [1.25794774e+02, -5.83537876e-04, -2.57431876e-03,  6.19251463e-01,
                                               1.38783126e+02,  1.03390335e-01, -2.12537252e-02, -5.28821722e-01,
                                               1.29719019e+02, -4.56132270e-04])

assert np.allclose(refractive.flatten()[::21], [1.25857198e+02, -5.89476000e-04, -2.56123378e-03,  6.19115458e-01,
                                                1.38850488e+02,  1.03201081e-01, -2.10969520e-02, -5.31725056e-01,
                                                1.29760585e+02, -4.48233446e-04])

assert not np.allclose(geometric, refractive)
