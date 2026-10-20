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
            [2.276512797010e+02, 4.247444534542e-04, 1.106734107879e-04, 5.609136825197e-02],
            [2.307288622021e+02, 6.580991399971e-04, 1.716517753891e-04, 6.943198305558e-02],
            [2.346717047728e+02, 1.164984787181e-03, 3.043583152191e-04, 9.210734119502e-02],
            [2.402265080250e+02, 2.607868940089e-03, 6.834658948378e-04, 1.381001038287e-01],
            [2.496499312861e+02, 9.723335071968e-03, 2.572463354838e-03, 2.662520389182e-01],
            [2.098472421771e+02, 2.413392785701e+01, 2.234597452017e+00, 5.387278008584e-06],
            [2.496493645983e+02, 9.726865161284e-03, 2.573426326947e-03, -2.663284227327e-01],
            [2.402254168113e+02, 2.609668569079e-03, 6.839448255818e-04, -1.381756012516e-01],
            [2.346701299283e+02, 1.166181685193e-03, 3.046742194925e-04, -9.218226792612e-02],
            [2.307268146713e+02, 6.590009454186e-04, 1.718888005229e-04, -6.950688130959e-02],
            [2.276487692701e+02, 4.254726536743e-04, 1.108643266718e-04, -5.616615269498e-02],
        ]
    ),
), "Magnus Zeeman spectral radiance has drifted from the expected values"
