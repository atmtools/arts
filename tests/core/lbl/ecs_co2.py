import pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.Workspace()

ws.abs_speciesSet(species=["CO2-626"])

ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")

ws.absorption_bandsFromAbsorbtionLines()

ws.jacobian_targets = pyarts.arts.JacobianTargets()
ws.select_abs_species = []  # All species
ws.atm_pointInit()
ws.atm_point.temperature = 295  # At room temperature
ws.atm_point.pressure = 1e5
ws.atm_point[ws.abs_species[0]] = 0.21  # At 21% Oxygen
ws.atm_point[pyarts.arts.SpeciesEnum("N2")] = 0.79  # At 21% Oxygen
ws.atm_point.mag = [40e-6, 20e-6, 10e-6]

ws.jacobian_targetsInit()
ws.Wigner6Init()

ws.ecs_dataInitNEWNEW()
ws.ecs_dataAddTran2011NEWNEW()
ws.ecs_dataAddRodrigues1997NEWNEW()
ws.ecs_dataAddMeanAirNEWNEW(vmrs=[0.21, 0.79], species=["O2", "N2"])

x = pyarts.arts.ArrayOfIndex()
ws.SortedQuantumIdentifiersOfBands(
    x, criteria="IntegratedIntensity", reverse=True
)

f2c = pyarts.arts.convert.freq2kaycm

y = pyarts.arts.AbsorptionBands(ws.absorption_bands)

for i in range(len(x)):
    ws.absorption_bands = [y[x[i]]]

    ws.f_grid = np.linspace(
        ws.absorption_bands[0].data.lines[0].f0 * 0.9,
        ws.absorption_bands[0].data.lines[-1].f0 * 1.1,
        10001,
    )  # around the band

    plt.figure(1, figsize=(10, 10))
    plt.clf()
    ws.absorption_bands[0].data.lineshape = "VP_ECS_HARTMANN"
    ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
    ws.propmat_clearskyAddLines2()
    plt.semilogy(f2c(ws.f_grid), ws.propmat_clearsky[:].T[0])

    ws.absorption_bands[0].data.lineshape = "VP_LTE"
    ws.propmat_clearskyInit(propmat_clearsky_agenda_checked=1)
    ws.propmat_clearskyAddLines2()
    plt.semilogy(f2c(ws.f_grid), ws.propmat_clearsky[:].T[0])
    plt.title(f"band {i}; qid: {ws.absorption_bands[0].key}")
    plt.ylabel("Absorption [1/m]")
    plt.xlabel("Frequency [cm-1]")
    plt.legend(["Line mixing ON", "Line mixing OFF"])
    plt.show()
