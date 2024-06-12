import pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.Workspace()

ws.absorption_speciesSet(species=["O2-66", "H2O-PWR98"])
ws.ReadCatalogData()

bandkey = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"
ws.absorption_bandsSelectFrequency(fmax=120e9)
ws.absorption_bandsKeepID(id=bandkey)

t = []
for a in ws.absorption_bands[0].data.lines:
    if a.f0 > 5e9 and a.f0 < 120e9:
        t.append(a)
ws.absorption_bands[0].data.lines = t

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

ws.absorption_bands[0].data.lineshape = "VP_ECS_MAKAROV"
ws.propagation_matrixInit()
ws.propagation_matrixAddLines()
ws.propagation_matrixAddPredefined()

plt.clf()
plt.semilogy(ws.frequency_grid/1e9, ws.propagation_matrix[:, 0])

ws.ReadCatalogData()

bandkey = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"
ws.absorption_bandsSelectFrequency(fmax=1200e9)
ws.absorption_bandsKeepID(id=bandkey)

t = []
for a in ws.absorption_bands[0].data.lines:
    if a.f0 > 5e9 and a.f0 < 120e9:
        t.append(a)
ws.absorption_bands[0].data.lines = t

ws.propagation_matrixInit()
ws.propagation_matrixAddLines()
ws.propagation_matrixAddPredefined()

plt.semilogy(ws.frequency_grid/1e9, ws.propagation_matrix[:, 0], ':')

for line in ws.absorption_bands[0].data.lines:
    line.ls.remove("Y")
    line.ls.remove("G")
    line.ls.remove("DV")

ws.propagation_matrixInit()
ws.propagation_matrixAddLines()
ws.propagation_matrixAddPredefined()

plt.semilogy(ws.frequency_grid/1e9, ws.propagation_matrix[:, 0])


ws.absorption_speciesSet(species=["O2-PWR98", "H2O-PWR98"])
ws.ReadCatalogData()
ws.propagation_matrixInit()
ws.propagation_matrixAddPredefined()
plt.semilogy(ws.frequency_grid/1e9, ws.propagation_matrix[:, 0], ":")



plt.legend(["ECS", "ONLINE", "NO-LM", "PWR"])
