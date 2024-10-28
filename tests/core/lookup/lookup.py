import pyarts
import numpy as np
import matplotlib.pyplot as plt
from time import time

# %% Setup workspace

ws = pyarts.Workspace()

ws.absorption_speciesSet(species=["CO2-626"])

ws.ReadCatalogData()
for key in ws.absorption_bands:
    ws.absorption_bands[key].cutoff = "ByLine"
    ws.absorption_bands[key].cutoff_value = 750e9

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field["t"] = 295.0

ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)

v = np.linspace(400, 2500, 1001)
ws.frequency_grid = pyarts.arts.convert.kaycm2freq(v)

t = time()
ws.absorption_lookup_tableCalc(
    longitude=0, latitude=0, temperature_perturbation=np.linspace(-20, 20, 11)
)
print (time() - t, "s to train the LUT")

ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")

pos = [100e3, 0, 0]
los = [120.0, 0.0]
ws.ray_pathGeometric(pos=pos, los=los, max_step=500.0)

# %% Checks and settings for LBL

t = time()
ws.propagation_matrix_agendaAuto(use_absorption_lookup_table=0)
ws.spectral_radianceClearskyEmission()
lbl = ws.spectral_radiance[:, 0] * 1.0
print (time() - t, "s to compute the LBL spectral radiance")

# %% Checks and for lookup calculations

t = time()
ws.propagation_matrix_agendaAuto(use_absorption_lookup_table=1)
ws.spectral_radianceClearskyEmission()
lut = ws.spectral_radiance[:, 0] * 1.0
print (time() - t, "s to compute the LUT spectral radiance")

# %% Plot the results

plt.semilogy(v, lbl, label="LBL")
plt.semilogy(v, lut, ":", label="LUT")
plt.legend()
plt.xlabel("Frequency [cm$^{-1}$]")
plt.ylabel("Radiance [W/m$^2$/sr/Hz]")

# %% Check the difference

assert np.allclose(lut, lut)
