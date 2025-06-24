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
ws.absorption_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldPlanet(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Add a sun

ws.sunBlackbody()
ws.sunsAddSun(suns=ws.suns)

# %% Checks and settings

ws.spectral_radiance_transform_operatorSet(option="Tb")
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
            [5.33326282e03, 5.33326282e03, 5.33326282e03],
            [2.73804755e00, 2.73801974e00, 2.73802787e00],
            [2.74724222e00, 2.74702051e00, 2.74708569e00],
            [2.76264051e00, 2.76189614e00, 2.76211620e00],
            [2.78427078e00, 2.78251808e00, 2.78303923e00],
            [2.81213495e00, 2.80873894e00, 2.80975463e00],
            [2.84621135e00, 2.84039654e00, 2.84214599e00],
        ]
    ),
)
