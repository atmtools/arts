import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

# Download catalogs
pyarts.data.download()

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range
line_f0 = 118750348044.712
ws.freq_grid = [line_f0]

# %% Species and line absorption
ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.abs_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.spectral_propmat_agendaAuto()

# %% Grids and planet
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atm_fieldIGRF(time="2000-03-11 14:39:37")

# %% Add a sun
ws.sunBlackbody()
ws.suns = [ws.sun]

# %% Checks and settings
ws.spectral_rad_transform_operatorSet(option="Tb")
ws.spectral_rad_space_agendaSet(option="SunOrCosmicBackground")
ws.spectral_rad_surface_agendaSet(option="Blackbody")

# %% Core calculations
pos = [90e3, 0, 0]
zens = np.linspace(0, 2, 21)
azis = np.linspace(-180, 180, 21)
res = np.empty((len(zens), len(azis)))
for izen in range(len(zens)):
    for iazi in range(len(azis)):
        los = [zens[izen], azis[iazi]]
        ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)
        ws.spectral_radClearskyEmission()
        ws.spectral_radApplyUnitFromSpectralRadiance()
        res[izen, iazi] = ws.spectral_rad[0][0]

# FIXME: Use some sort of Imager for measurement_vec for the above

r, theta = np.meshgrid(zens, np.rad2deg(azis))
fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
ax.contourf(theta, r, res.T)

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

assert np.allclose(
    res[::3, ::7],
    np.array(
        [[5333.26510189, 5333.26510189, 5333.26510189],
         [17.54346643,   17.42489325,   17.45953778],
         [17.61367450,   17.37693775,   17.44655483],
         [17.68367603,   17.32918821,   17.43410485],
         [17.75346986,   17.28164635,   17.42218858],
         [17.82305482,   17.23431393,   17.41080683],
         [17.89242977,   17.18719274,   17.39996041]]
    ),
)
