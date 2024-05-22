import pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

ws.frequency_grid = [pyarts.arts.convert.wavelen2freq(700e-9)]

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
ws.ray_path_observer_agendaSet(option="Geometric")

# %% Core calculations

pos = [90e3, 0, 0]
zas = np.linspace(0, 5, 21)
aas = np.linspace(-180, 180, 21)
res = np.empty((len(zas), len(aas)))
for iza in range(len(zas)):
    for iaa in range(len(aas)):
        los = [zas[iza], aas[iaa]]
        ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)
        ws.ray_path_suns_pathFromPathObserver()
        ws.spectral_radianceClearskyRayleighScattering()
        ws.spectral_radianceApplyUnitFromSpectralRadiance()
        res[iza, iaa] = ws.spectral_radiance[0][0]

# FIXME: Use some sort of Imager for measurement_vector for the above

r, theta = np.meshgrid(zas, np.rad2deg(aas))
fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
ax.contourf(theta, r, res.T)

pos = [90e3, 0, 0]
za = 18.0
aa = -90.0
ws.ray_pathGeometric(pos=pos, los=[za, aa], max_step=1000.0)
ws.ray_path_suns_pathFromPathObserver()
ws.spectral_radianceClearskyRayleighScattering()
ws.spectral_radianceApplyUnitFromSpectralRadiance()


assert np.allclose(
    res.flatten()[::10],
    np.array(
        [
            5772.0,
            5772.0,
            5772.0,
            5772.0,
            5772.0,
            644.09429354,
            644.09429354,
            644.09429354,
            644.09429354,
            644.09429354,
            644.09429354,
            644.09429354,
            644.09429354,
            644.09429354,
            644.09429354,
            644.09429354,
            644.09429354,
            644.08595839,
            644.08595839,
            644.08595839,
            644.08595839,
            644.08595839,
            644.08595839,
            644.08595839,
            644.10262547,
            644.10262547,
            644.10262547,
            644.10262547,
            644.10262547,
            644.10262547,
            644.11927965,
            644.11927965,
            644.08595839,
            644.08595839,
            644.06093357,
            644.06093357,
            644.03587962,
            644.03587962,
            644.06927841,
            644.06927841,
            644.03587962,
            644.03587962,
            644.06927841,
            644.06927841,
            644.06927841,
        ]
    ),
)
