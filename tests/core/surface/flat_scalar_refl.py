import pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

ws.frequency_grid = np.linspace(22e9, 23e9, 1001)

# %% Species and line absorption

ws.absorption_speciesSet(species=["H2O-PWR98", "O2-PWR98"])
ws.ReadCatalogData()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 300.0
ws.surface_field["flat scalar reflectance"] = 0.99
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)

# %% Checks and settings

ws.spectral_radiance_observer_agendaSet(option="EmissionNoSensor")
ws.spectral_radiance_surface_agendaSet(option="FlatScalarReflectance")
ws.ray_path_observer_agendaSetGeometric()
ws.spectral_radiance_transform_operatorSet(option="RJBT")

# %% Position

pos = [0, 0, 0]

# %% Up calculations

los = [0.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)
ws.spectral_radianceClearskyEmission()
ws.spectral_radianceApplyUnitFromSpectralRadiance()
up = ws.spectral_radiance[:, 0] * 1.0

# %% Down calculations

los = [180.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)
ws.spectral_radianceClearskyEmission()
ws.spectral_radianceApplyUnitFromSpectralRadiance()
down = ws.spectral_radiance[:, 0] * 1.0

# %% Show results

plt.plot((ws.frequency_grid/1e9), up)
plt.plot((ws.frequency_grid/1e9), down)
plt.xlabel("Frequency offset [MHz]")
plt.ylabel("Spectral radiance [K]")

# %% Check but only the basics

assert np.all(up < down)
assert np.all(down - up < 3.0)
assert np.allclose((down-up)[::100], [2.29985014, 2.29092721, 2.28374728, 2.2788642, 2.27601317,
                                      2.27506843, 2.27599661, 2.27867851, 2.28296023, 2.28867035,
                                      2.29563091])
