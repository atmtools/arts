import pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.frequency_grid = [line_f0]

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequency(fmin=40e9, fmax=120e9, by_line=1)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Add a sun

ws.sunBlackbody()
ws.sunsAddSun(suns=ws.suns)

# %% Checks and settings

ws.spectral_radiance_unit = "Tb"
ws.spectral_radiance_space_agendaSet(option="SunOrCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

# %% Core calculations

pos = [90e3, 0, 0]
zas = np.linspace(0, 2, 21)
aas = np.linspace(-180, 180, 21)
res = np.empty((len(zas), len(aas)))
for iza in range(len(zas)):
    for iaa in range(len(aas)):
        los = [zas[iza], aas[iaa]]
        ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)
        ws.spectral_radianceClearskyEmission()
        ws.spectral_radianceApplyUnitFromSpectralRadiance()
        res[iza, iaa] = ws.spectral_radiance[0][0]

# FIXME: Use some sort of Imager for measurement_vector for the above

r, theta = np.meshgrid(zas, np.rad2deg(aas))
fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
ax.contourf(theta, r, res.T)

assert np.allclose(
    res[::3, ::7],
    np.array(
        [
            [5749.74670977, 5673.9300261, 5420.98020555],
            [7.13171303, 8.7751568, 14.49545782],
            [12.19459019, 12.87760697, 15.50152675],
            [17.06697185, 17.10784627, 17.28252788],
            [20.24466269, 20.08448805, 19.18630347],
            [20.87660536, 20.95869027, 20.56264444],
            [18.80695427, 19.50096608, 20.97105106],
        ]
    ),
)
