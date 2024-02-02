import pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

# %% Sensor

ws.frequency_grid = np.linspace(-50e6, 50e6, 101) + 118750348044.712

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66-118e9-119e9"])
ws.abs_lines_per_species = pyarts.arts.ArrayOfArrayOfAbsorptionLines()
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(
    ws.abs_lines_per_species, basename="lines/"
)
ws.absorption_bandsFromAbsorbtionLines(
    abs_lines_per_species=ws.abs_lines_per_species
)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.options.SurfaceKey("t")] = 295.0
t = pyarts.arts.GriddedField3.fromxml("planets/Earth/afgl/tropical/t.xml")
ws.atmospheric_fieldInit(toa=100e3)
ws.atmospheric_fieldAddGriddedData(
    key=pyarts.arts.String("t"),
    data=pyarts.arts.GriddedField3.fromxml(
        "planets/Earth/afgl/tropical/t.xml"
    ),
)
ws.atmospheric_fieldAddGriddedData(
    key=pyarts.arts.String("p"),
    data=pyarts.arts.GriddedField3.fromxml(
        "planets/Earth/afgl/tropical/p.xml"
    ),
)
ws.atmospheric_field[pyarts.arts.SpeciesEnum("O2")] = 0.21
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_space_agendaSet()
ws.spectral_radiance_surface_agendaSet()

# %% Core calculations
ws.propagation_pathGeometric(pos=[50e3, 0, 0], los=[0, 0], max_step=1e2)
ws.spectral_radianceStandardEmission()


import matplotlib.pyplot as plt
alts = np.linspace(0, 1e5, 1001)
s = pyarts.arts.SpectralRadianceOperator()
ws.spectral_radiance_operatorGeometricPlanar(s, altitude_grid=alts)
srad = s.geometric_planar(ws.frequency_grid, [50e3,0,0], [180,180])
plt.semilogy(ws.frequency_grid, ws.spectral_radiance)
plt.show()
plt.semilogy(ws.frequency_grid, srad)
