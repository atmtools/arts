import pyarts
import numpy as np

ws = pyarts.Workspace()

ws.absorption_speciesSet(species=["O2-66"])

ws.abs_lines_per_species = pyarts.arts.ArrayOfArrayOfAbsorptionLines()
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(
    ws.abs_lines_per_species, basename="lines/"
)

bandkey = "O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"

ws.absorption_bandsFromAbsorbtionLines(
    abs_lines_per_species=ws.abs_lines_per_species
)
ws.absorption_bandsSelectFrequency(fmax=120e9)
ws.absorption_bandsKeepID(id=bandkey)

t = []
for a in ws.absorption_bands[0].data.lines:
    if a.f0 > 5e9 and a.f0 < 120e9:
        t.append(a)
ws.absorption_bands[0].data.lines = t

ws.WignerInit()
ws.frequency_grid = np.linspace(40e9, 130e9, 10001)  # around the band

ws.jacobian_targets = pyarts.arts.JacobianTargets()
ws.atmospheric_pointInit()
ws.atmospheric_point.temperature = 295  # At room temperature
ws.atmospheric_point.pressure = 1e5
ws.atmospheric_point[pyarts.arts.SpeciesEnum("O2")] = 0.21  # At 21% Oxygen
ws.atmospheric_point[pyarts.arts.SpeciesEnum("N2")] = 0.79  # At 79% Nitrogen
ws.atmospheric_point.mag = [40e-6, 20e-6, 10e-6]
ws.propagation_path_point

ws.jacobian_targetsInit()

ws.ecs_dataInit()
ws.ecs_dataAddMakarov2020()
ws.ecs_dataAddMeanAir(vmrs=[1], species=["N2"])

ws.absorption_bands[0].data.lineshape = "VP_ECS_MAKAROV"
ws.propagation_matrixInit()
ws.propagation_matrixAddLines()
