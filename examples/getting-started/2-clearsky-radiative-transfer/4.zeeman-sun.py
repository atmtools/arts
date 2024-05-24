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
    res.flatten()[::10],
    np.array(
        [
            5747.62549888,
            5247.40532663,
            5747.62549888,
            5289.73738873,
            5562.37239443,
            5576.07412529,
            5282.58761735,
            7.48965076,
            20.96754703,
            14.72144364,
            20.94041283,
            19.57450237,
            20.25807923,
            20.51723219,
            16.78579824,
            19.83300643,
            14.47034222,
            18.03924429,
            15.81215198,
            17.3353868,
            17.43013079,
            18.44360588,
            16.88784963,
            18.44360588,
            17.8380293,
            19.1207752,
            19.95138437,
            17.80859736,
            20.76520606,
            13.54010971,
            20.25560369,
            9.86011333,
            16.13606293,
            13.89721384,
            7.88371635,
            20.2061116,
            7.51673719,
            20.21309617,
            17.44397868,
            19.24973875,
            20.63782912,
            20.97577602,
            16.07983693,
            11.21084166,
            16.07983693,
        ]
    ),
)
