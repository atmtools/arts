import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

# Download catalogs
pyarts.data.download()

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range
ws.frequency_grid = [pyarts.arts.convert.wavelen2freq(700e-9)]

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
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Add a sun
ws.sunBlackbody()
ws.suns = [ws.sun]

# %% Checks and settings
ws.spectral_radiance_transform_operator = "Tb"
ws.spectral_radiance_space_agendaSet(option="SunOrCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.ray_path_observer_agendaSetGeometric()
ws.propagation_matrix_scattering_agendaSet(option="AirSimple")

# %% Core calculations
pos = [90e3, 0, 0]
zas = np.linspace(0, 5, 21)
aas = np.linspace(-180, 180, 21)
res = np.empty((len(zas), len(aas)))
for iza in range(len(zas)):
    for iaa in range(len(aas)):
        los = [zas[iza], aas[iaa]]
        ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)
        ws.ray_path_suns_pathFromPathObserver(just_hit=1)
        ws.spectral_radianceClearskyRayleighScattering()
        ws.spectral_radianceApplyUnitFromSpectralRadiance()
        res[iza, iaa] = ws.spectral_radiance[0][0]

# FIXME: Use some sort of Imager for measurement_vector for the above

r, theta = np.meshgrid(zas, np.rad2deg(aas))
fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
ax.contourf(theta, r, res.T)

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

assert np.allclose(
    res[1::3, 1::7],
    np.array([[5771.99991010, 5771.99991010, 5771.99991010],
              [ 643.69008277,  643.69008277,  643.69008277],
              [ 643.69030757,  643.69030764,  643.69030764],
              [ 643.69065708,  643.69065715,  643.69065728],
              [ 643.69113428,  643.69113428,  643.69113441],
              [ 643.69174257,  643.69174270,  643.69174290],
              [ 643.69248695,  643.69248708,  643.69248747]]),
)
