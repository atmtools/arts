import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.Workspace()

ws.abs_speciesSet(species=["O2-66", "H2O-PWR98"])
ws.ReadCatalogData()

bandkey = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"
ws.abs_bandsSelectFrequencyByBand(fmax=120e9)
ws.abs_bandsKeepID(id=bandkey)

t = []
for a in ws.abs_bands[bandkey].lines:
    if a.f0 > 5e9 and a.f0 < 120e9:
        t.append(a)
ws.abs_bands[bandkey].lines = t


def calc(ws, lineshape=None):
    if lineshape is not None:
        ws.abs_bands[bandkey].lineshape = lineshape
    ws.propagation_matrixInit()
    ws.propagation_matrixAddLines()
    ws.propagation_matrixAddPredefined()
    return 1.0 * ws.propagation_matrix[:, 0]


ws.WignerInit()
ws.frequency_grid = np.linspace(20e9, 140e9, 10001)  # around the band

ws.jacobian_targets = pyarts.arts.JacobianTargets()
ws.atmospheric_pointInit()
ws.atmospheric_point.temperature = 295  # At room temperature
ws.atmospheric_point.pressure = 1e5
ws.atmospheric_point[pyarts.arts.SpeciesEnum("CO2")] = 0.0
ws.atmospheric_point[pyarts.arts.SpeciesEnum("liquidcloud")] = 0.0
ws.atmospheric_point[pyarts.arts.SpeciesEnum("O2")] = 0.21  # At 21% Oxygen
ws.atmospheric_point[pyarts.arts.SpeciesEnum("H2O")] = 0.01  # At 1% Water
ws.atmospheric_point[pyarts.arts.SpeciesEnum("N2")] = 0.79  # At 79% Nitrogen
ws.atmospheric_point.mag = [40e-6, 20e-6, 10e-6]
ws.ray_path_point

ws.jacobian_targetsInit()

ws.ecs_dataInit()
ws.ecs_dataAddMakarov2020()
ws.ecs_dataAddMeanAir(vmrs=[1], species=["N2"])


plt.clf()

# Online data
plt.semilogy(ws.frequency_grid / 1e9, calc(ws), label="Online", lw=3)

# ECS data
plt.semilogy(ws.frequency_grid / 1e9, calc(ws, "VP_ECS_MAKAROV"), label="ECS")

# Remove line mixing
ws.abs_bands.clear_linemixing()
plt.semilogy(
    ws.frequency_grid / 1e9, calc(ws, "VP_LTE"), label="No linemixing"
)

# 1st order line mixing
ws.abs_bands[bandkey].lineshape = "VP_ECS_MAKAROV"
ws.abs_bandsLineMixingAdaptation(
    temperatures=np.linspace(200, 350, 16),
    band_key=bandkey,
    rosenkranz_fit_order=1,
)
plt.semilogy(
    ws.frequency_grid / 1e9,
    calc(ws, "VP_LTE"),
    "--",
    label="1st Order Rosenkranz",
)

# 2nd order line mixing
ws.abs_bands[bandkey].lineshape = "VP_ECS_MAKAROV"
ws.abs_bandsLineMixingAdaptation(
    temperatures=np.linspace(200, 350, 16),
    band_key=bandkey,
    rosenkranz_fit_order=2,
)
plt.semilogy(
    ws.frequency_grid / 1e9,
    calc(ws, "VP_LTE"),
    ":",
    label="2nd Order Rosenkranz",
)

# Using PWR98
ws.abs_speciesSet(species=["O2-PWR98", "H2O-PWR98"])
ws.ReadCatalogData()
plt.semilogy(ws.frequency_grid / 1e9, calc(ws), label="PWR98")

plt.legend()

## Test here if we can save and load the workspace
ws.savexml("test.xml")

ws2 = pyarts.Workspace.fromxml("test.xml")
assert np.all(ws.atmospheric_point.mag == ws2.atmospheric_point.mag)
