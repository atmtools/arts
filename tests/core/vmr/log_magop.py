import os

import pyarts3 as pyarts
import numpy as np
from copy import copy

NF = 1001
noise = 1

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

ws.freq_grid = np.linspace(10e9, 400e9, NF)

# %% Species and line absorption

ws.abs_speciesSet(species=["H2O-PWR98", "O2-PWR98"])
ws.ReadCatalogData()
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.spectral_propmat_agendaAuto()

# %% Grids and planet

ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)

# %% Checks and settings

ws.spectral_rad_transform_operatorSet(option="Tb")
ws.ray_path_observer_agendaSetGeometric()
ws.rte_option = "magop"

pos = [0e3, 0, 0]
los = [20.0, 0.0]
ws.measurement_sensorSimple(pos=pos, los=los)

RAT = 0.8
fieldg = copy(ws.atm_field["H2O"])

fieldg.data /= RAT

ws.RetrievalInit()
ws.RetrievalAddSpeciesVMR(species="H2O", matrix=np.diag(np.ones((50)) * 1e-4))
ws.RetrievalFinalizeDiagonal()

ws.jac_targetsToggleLogarithmicAtmTarget(key="H2O")

ws.measurement_vecFromSensor()
transmat = np.array(ws.measurement_vec, copy=True)

ws.atm_field["H2O"] = fieldg
ws.model_state_vec_aprioriFromData()
ws.measurement_vecFromSensor()
apriori = np.array(ws.measurement_vec, copy=True)

ws.measurement_vec_error_covmatConstant(value=noise**2)
ws.measurement_vec = transmat
ws.OEM(method="lm", lm_ga_settings=[10, 2, 2, 100, 1, 99])
vmr_jacobian = np.array(ws.measurement_vec_fit, copy=True)

np.testing.assert_array_less(
    np.linalg.norm(transmat - vmr_jacobian),
    np.linalg.norm(transmat - apriori),
)

if "ARTS_HEADLESS" not in os.environ:
    import matplotlib.pyplot as plt

    plt.plot(ws.freq_grid / 1e9, transmat, label="transmission target")
    plt.plot(ws.freq_grid / 1e9, apriori, label="apriori")
    plt.plot(ws.freq_grid / 1e9, vmr_jacobian, label="VMR-Jacobian fit")
    plt.xlabel("Frequency [GHz]")
    plt.ylabel("Brightness temperature [K]")
    plt.legend()
    plt.show()
