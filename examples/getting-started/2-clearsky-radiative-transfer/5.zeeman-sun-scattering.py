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
ws.propagation_matrix_scattering_agendaSet(option="AirSimple")

# %% Core calculations
pos = [90e3, 0, 0]
zas = np.linspace(0, 5, 21)
aas = np.linspace(-180, 180, 21)
res = np.empty((len(zas), len(aas)))
for iza in range(len(zas)):
    for iaa in range(len(aas)):
        los = [zas[iza], aas[iaa]]
        ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)
        ws.ray_path_suns_pathFromPathObserver(just_hit=1)
        ws.spectral_radianceClearskyRayleighScattering()
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
            5771.9999101,
            5771.9999101,
            5771.9999101,
            5771.9999101,
            5771.9999101,
            643.69000049,
            643.69000049,
            643.69003444,
            643.69003444,
            643.69008219,
            643.69008219,
            643.69014321,
            643.69014321,
            643.69021817,
            643.69021817,
            643.69030699,
            643.69030699,
            643.69040942,
            643.69040942,
            643.69052583,
            643.69052583,
            643.69065631,
            643.69065631,
            643.69065631,
            643.69080123,
            643.69080123,
            643.69096008,
            643.69096008,
            643.69113369,
            643.69113369,
            643.69132195,
            643.69132195,
            643.69152465,
            643.69152465,
            643.69174225,
            643.69174225,
            643.69197481,
            643.69197481,
            643.69222278,
            643.69222278,
            643.69248637,
            643.69248637,
            643.69276562,
            643.69276562,
            643.69276562,
        ]
    ),
)
