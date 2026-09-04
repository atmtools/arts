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
ws.rte_option = "magop"
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
            [2.276512821092e02, 4.261021957888e-04, 1.027468021420e-04, 5.687927378416e-02],
            [2.307288660365e02, 6.602011364976e-04, 1.593711565846e-04, 7.040747139975e-02],
            [2.346717117551e02, 1.168703614880e-03, 2.826220439971e-04, 9.340173537775e-02],
            [2.402265237320e02, 2.616196566095e-03, 6.348256699340e-04, 1.400415316211e-01],
            [2.496499814427e02, 9.754569416941e-03, 2.391349108080e-03, 2.699971405358e-01],
            [2.099028163397e02, 2.416037218539e01, 1.740559369663e00, 5.523843341848e-06],
            [2.496494147479e02, 9.758111279439e-03, 2.392246667185e-03, -2.700745982318e-01],
            [2.402254325160e02, 2.618002063414e-03, 6.352711243096e-04, -1.401180900068e-01],
            [2.346701369088e02, 1.169904391432e-03, 2.829156569959e-04, -9.347771470715e-02],
            [2.307268185042e02, 6.611058571622e-04, 1.595913768142e-04, -7.048342160883e-02],
            [2.276487716773e02, 4.268327476638e-04, 1.029241438602e-04, -5.695510840674e-02],
        ]
    ),
), "Magnus Zeeman spectral radiance has drifted from the expected values"
