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
            [2.276551235399e02, 4.259411802252e-04, 1.027078921965e-04, 5.686609432459e-02],
            [2.307318484287e02, 6.600222918109e-04, 1.593279395564e-04, 7.038684472572e-02],
            [2.346737153328e02, 1.168587845157e-03, 2.825940092634e-04, 9.338410284805e-02],
            [2.402280899194e02, 2.616549259585e-03, 6.349109213630e-04, 1.400459368601e-01],
            [2.496539453642e02, 9.761561980270e-03, 2.393061571382e-03, 2.701384329393e-01],
            [2.098955622505e02, 2.417326246041e01, 1.741633379836e00, 5.526631596808e-06],
            [2.496533784517e02, 9.765105367691e-03, 2.393959518344e-03, -2.702159091422e-01],
            [2.402269987255e02, 2.618354729037e-03, 6.353563693580e-04, -1.401224878107e-01],
            [2.346721406968e02, 1.169788358766e-03, 2.828875578587e-04, -9.346007246592e-02],
            [2.307298009667e02, 6.609267549322e-04, 1.595480962493e-04, -7.046280541539e-02],
            [2.276526123549e02, 4.266716225771e-04, 1.028852062035e-04, -5.694198344678e-02],
        ]
    ),
), "Linear-source Magnus Zeeman spectral radiance has drifted from the expected values"
