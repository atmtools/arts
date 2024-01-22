import pyarts
import numpy as np

ws = pyarts.Workspace()

ws.abs_speciesSet(species=["O2-66"])

ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")

bandkey = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"

ws.absorption_bandsFromAbsorbtionLines()
ws.absorption_bandsSelectFrequency(fmax=120e9)
ws.absorption_bandsKeepID(id=bandkey)

t = []
for a in ws.absorption_bands[0].data.lines:
    if a.f0 > 5e9 and a.f0 < 120e9:
        t.append(a)
ws.absorption_bands[0].data.lines = t

ws.Wigner6Init()
ws.f_grid = np.linspace(40e9, 130e9, 10001)  # around the band

ws.jacobian_targets = pyarts.arts.JacobianTargets()
ws.select_abs_species = []  # All species
ws.atm_pointInit()
ws.atm_point.temperature = 295  # At room temperature
ws.atm_point.pressure = 1e5
ws.atm_point[ws.abs_species[0]] = 0.21  # At 21% Oxygen
ws.atm_point[pyarts.arts.SpeciesEnum("N2")] = 0.79  # At 21% Oxygen
ws.atm_point.mag = [40e-6, 20e-6, 10e-6]

ws.jacobian_targetsInit()

ws.ecs_dataInitNEWNEW()
ws.ecs_dataAddMakarov2020NEWNEW()
ws.ecs_dataAddMeanAirNEWNEW(vmrs=[1], species=["N2"])

ws.absorption_bands[0].data.lineshape = "VP_ECS_MAKAROV"
ws.propmat_clearskyInit()
ws.propmat_clearskyAddLines2()
