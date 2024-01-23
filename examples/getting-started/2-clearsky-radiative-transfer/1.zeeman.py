import pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

# %% Sensor

ws.f_grid = np.linspace(-4e6, 4e6, 1001) + 118750348044.712

# %% Species and line absorption

ws.abs_speciesSet(species=["O2-Z-66-118e9-120e9"])
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")
ws.Wigner6Init()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.options.SurfaceKey("t")] = 295.0
t = pyarts.arts.GriddedField3.fromxml(
    "planets/Earth/afgl/tropical/t.xml"
)
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
ws.atmospheric_field[ws.abs_species[0]] = 0.21
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_background_space_agendaSet()
ws.spectral_radiance_background_surface_agendaSet()

# %% Core calculations
ws.propagation_pathGeometric(pos=[300e3, 0, 0], los=[180, 0])
ws.spectral_radianceStandardEmission()
