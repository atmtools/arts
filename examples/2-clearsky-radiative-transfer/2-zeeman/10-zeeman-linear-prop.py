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
ws.rte_option = "linprop"
ws.spectral_radClearskyEmission()
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
        [[ 2.27655095e+02,  4.24584408e-04,  1.10632212e-04,  5.60784161e-02],
         [ 2.30731819e+02,  6.57921523e-04,  1.71606497e-04,  6.94117038e-02],
         [ 2.34673684e+02,  1.16487040e-03,  3.04330863e-04,  9.20900445e-02],
         [ 2.40228051e+02,  2.60822199e-03,  6.83565997e-04,  1.38104591e-01],
         [ 2.49653876e+02,  9.73030171e-03,  2.57436351e-03,  2.66391630e-01],
         [ 2.09838939e+02,  2.41459106e+01,  2.23292127e+00,  5.39076306e-06],
         [ 2.49653309e+02,  9.73383330e-03,  2.57532696e-03, -2.66468032e-01],
         [ 2.40226960e+02,  2.61002159e-03,  6.84044937e-04, -1.38180081e-01],
         [ 2.34672109e+02,  1.16606703e-03,  3.04646705e-04, -9.21649616e-02],
         [ 2.30729772e+02,  6.58823071e-04,  1.71843459e-04, -6.94866124e-02],
         [ 2.27652584e+02,  4.25312499e-04,  1.10823101e-04, -5.61532543e-02]]
    ),
), "Values have drifted from expected results in spectral radiance"
