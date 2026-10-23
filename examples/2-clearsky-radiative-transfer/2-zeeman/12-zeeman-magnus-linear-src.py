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
ws.rte_option = "magop_linsrc"
ws.spectral_radClearskyEmission()
ws.spectral_radApplyUnitFromSpectralRadiance()

# %% Show results
fig, ax = pyarts.plot(ws.spectral_rad, freqs=(ws.freq_grid - line_f0) / 1e6)
[a.set_xlabel("Frequency offset [MHz]") for a in ax.flatten()]
[a.set_ylabel("Spectral radiance [K]") for a in ax.flatten()]
fig.suptitle(f"Zeeman effect of {round(line_f0 / 1e6)} MHz O$_2$ line")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% Test

assert np.allclose(
    ws.spectral_rad[::100],
    np.array(
        [
            [2.276551211280e+02, 4.245839535741e-04, 1.106315156818e-04, 5.607837140657e-02],
            [2.307318445869e+02, 6.579208679653e-04, 1.716052439633e-04, 6.941164201663e-02],
            [2.346737083403e+02, 1.164869392326e-03, 3.043281406008e-04, 9.208995270599e-02],
            [2.402280742008e+02, 2.608220533943e-03, 6.835577220556e-04, 1.381044473892e-01],
            [2.496538951196e+02, 9.730305338849e-03, 2.574305681449e-03, 2.663913683550e-01],
            [2.098399613348e+02, 2.414680054663e+01, 2.235931391076e+00, 5.389997696170e-06],
            [2.496533282142e+02, 9.733836947902e-03, 2.575269070036e-03, -2.664677704185e-01],
            [2.402269830092e+02, 2.610020135222e-03, 6.840366455663e-04, -1.381799374798e-01],
            [2.346721337059e+02, 1.166066028219e-03, 3.046439756389e-04, -9.216486986011e-02],
            [2.307297971263e+02, 6.588224166393e-04, 1.718422009418e-04, -6.948655060651e-02],
            [2.276526099440e+02, 4.253120446975e-04, 1.108224017514e-04, -5.615320959551e-02],
        ]
    ),
), "Linear-source Magnus Zeeman spectral radiance has drifted from the expected values"
