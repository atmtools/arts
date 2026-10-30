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
ws.atm_fieldIGRF(time="2000-03-11 14:39:37")

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

assert np.allclose(geometric.flatten()[::21], [1.25485848e+02, -6.01567897e-04, -2.41577985e-03,  6.08670834e-01,
                                               1.38412669e+02,  1.10187275e-01, -2.02820836e-02, -5.14978283e-01,
                                               1.29364509e+02, -4.84799820e-04])

assert np.allclose(refractive.flatten()[::21], [1.25485844e+02, -6.01567599e-04, -2.41578138e-03,  6.08670881e-01,
                                                1.38412666e+02,  1.10187297e-01, -2.02820987e-02, -5.14978187e-01,
                                                1.29364505e+02, -4.84799527e-04])

# The line dispersion alone only changes the path slightly
assert 0 < np.abs(geometric - refractive).max() < 1e-4
