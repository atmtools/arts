import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt


def calc(ws, lineshape=None):
    if lineshape is not None:
        ws.abs_bands[bandkey].lineshape = lineshape
    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines()
    return 1.0 * ws.propagation_matrix[:, 0]


ws = pyarts.Workspace()
ws.abs_speciesSet(species=["CO2-626"])
ws.ReadCatalogData()

p = 1e5
ws.jac_targets = pyarts.arts.JacobianTargets()
ws.atm_pointInit()
ws.atm_point.temperature = 295  # At room temperature
ws.atm_point.pressure = p
ws.atm_point[pyarts.arts.SpeciesEnum("CO2")] = 400e-6
ws.atm_point[pyarts.arts.SpeciesEnum("O2")] = 0.21  # At 21% Oxygen
ws.atm_point[pyarts.arts.SpeciesEnum("N2")] = 0.79  # At 79% Nitrogen
ws.atm_point.mag = [40e-6, 20e-6, 10e-6]
ws.ray_path_point

ws.jac_targetsInit()
ws.WignerInit()

ws.ecs_dataInit()
ws.ecs_dataAddTran2011()
ws.ecs_dataAddRodrigues1997()
ws.ecs_dataAddMeanAir(vmrs=[0.21, 0.79], species=["O2", "N2"])

f2c = pyarts.arts.convert.freq2kaycm

bandkey = "CO2-626 ElecStateLabel X X kronigParity f f l2 2 2 parity NODEF + r 2 1 v1 1 0 v2 2 2 v3 1 0"
ws.abs_bands = {bandkey: ws.abs_bands[bandkey]}
ws.freq_grid = np.linspace(
    ws.abs_bands[bandkey].lines[0].f0 * 0.8,
    1.2 * ws.abs_bands[bandkey].lines[-1].f0,
    1001,
)

plt.clf()

# Online data
plt.semilogy(ws.freq_grid / 1e9, calc(ws), label="Online", lw=3)

# ECS data
plt.semilogy(ws.freq_grid / 1e9, calc(ws, "VP_ECS_HARTMANN"), label="ECS")

# Remove line mixing
ws.abs_bands.clear_linemixing()
plt.semilogy(
    ws.freq_grid / 1e9, calc(ws, "VP_LTE"), label="No linemixing"
)

# 1st order line mixing
ws.abs_bands[bandkey].lineshape = "VP_ECS_HARTMANN"
ws.abs_bandsLineMixingAdaptation(
    temperatures=np.linspace(200, 350, 16),
    band_key=bandkey,
    rosenkranz_fit_order=1,
)
plt.semilogy(
    ws.freq_grid / 1e9,
    calc(ws, "VP_LTE"),
    "--",
    label="1st Order Rosenkranz",
)

plt.legend()
