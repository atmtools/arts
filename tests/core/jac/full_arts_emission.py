import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.freq_grid = np.linspace(-20e6, 20e6, 101) + line_f0

# %% Species and line absorption

ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.abs_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atm_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_transform_operatorSet(option="Tb")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

grid = pyarts.arts.GriddedField3(
    name="VMR",
    data=np.ones((3, 1, 1)) * 0.2,
    grid_names=["Altitude", "Latitude", "Longitude"],
    grids=[[0, 50e3, 120e3], [0], [0]],
)

ws.atm_field[pyarts.arts.SpeciesEnum.O2] = grid

# %% Jacobian

ws.measurement_sensor = []
ws.jacobian_targetsInit()
ws.jacobian_targetsAddSpeciesVMR(species="O2")
ws.jacobian_targetsFinalize()

# %% Core calculations

pos = [100e3, 0, 0]
los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)

ws.spectral_radianceClearskyEmission()

x1 = 1.0 * np.array(ws.spectral_radiance)
dx1 = 1.0 * np.array(ws.spectral_radiance_jacobian)

DX = 1e-8
dx2 = []
for i in range(3):
    ws.atm_field[pyarts.arts.SpeciesEnum.O2].data = grid
    ws.atm_field[pyarts.arts.SpeciesEnum.O2].data.data[i] += DX

    ws.spectral_radianceClearskyEmission()

    x2 = 1.0 * np.array(ws.spectral_radiance)
    dx2.append((x2 - x1) / DX)
dx2 = np.array(dx2)

# Does this large discrepancy mean we need better code?
assert np.allclose(dx1 / dx2 - 1, 0, atol=0.02), "Should be within 2%"
