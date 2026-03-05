import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt
import time


def calc(ws, lineshape=None):
    time_start = time.time()
    if lineshape is not None:
        ws.abs_bands[bandkey].lineshape = lineshape
    ws.spectral_propmatInit()
    ws.spectral_propmatAddLines()
    print(f"Time to calculate spectral propmat: {round(1000*(time.time() - time_start), 3)} ms for lineshape {lineshape}")
    return 1.0 * ws.spectral_propmat[:, 0]

print("ECS CO2 tests")

ws = pyarts.Workspace()
ws.abs_speciesSet(species=["CO2-626"])

print("Reading catalog data...")
time_start = time.time()
ws.ReadCatalogData()
print(f"Time to read catalog data: {round(1000*(time.time() - time_start), 3)} ms")

p = 1e5
ws.jac_targets = pyarts.arts.JacobianTargets()
ws.atm_pointInit()
ws.atm_point.temperature = 295  # At room temperature
ws.atm_point.pressure = p
ws.atm_point[pyarts.arts.SpeciesEnum("CO2")] = 400e-6
ws.atm_point[pyarts.arts.SpeciesEnum("O2")] = 0.21  # At 21% Oxygen
ws.atm_point[pyarts.arts.SpeciesEnum("N2")] = 0.79  # At 79% Nitrogen
ws.atm_point.mag = [40e-6, 20e-6, 10e-6]
ws.ray_point

ws.jac_targetsInit()
ws.WignerInit()

ws.abs_ecs_dataInit()
ws.abs_ecs_dataAddTran2011()
ws.abs_ecs_dataAddRodrigues1997()
ws.abs_ecs_dataAddMeanAir(vmrs=[0.21, 0.79], species=["O2", "N2"])

f2c = pyarts.arts.convert.freq2kaycm

bandkey = "CO2-626 ElecStateLabel X X kronigParity f f r 2 1 l2 2 2 parity * + v1 4 0 v2 2 2 v3 1 0"
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

print ("Adaptation tests")

timestart = time.time()
# 1st order line mixing
ws.abs_bands[bandkey].lineshape = "VP_ECS_HARTMANN"
ws.abs_bandsLineMixingAdaptation(
    temperatures=np.linspace(200, 350, 16),
    band_key=bandkey,
    rosenkranz_fit_order=1,
)
print(f"Time for adaptation 1st order: {round(1000*(time.time() - timestart), 3)} ms")

plt.semilogy(
    ws.freq_grid / 1e9,
    calc(ws, "VP_LTE"),
    "--",
    label="1st Order Rosenkranz",
)

plt.legend()
