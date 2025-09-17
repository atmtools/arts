"""Spectral atmospheric flux operator"""

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts
from matplotlib import ticker

# %% Initialize the operator
#
fop = pyarts.recipe.SpectralAtmosphericFlux(
    species=["H2O-161", "O2-66", "N2-44", "CO2-626", "O3-XFIT"],
    remove_lines_percentile={"H2O": 70},
)

# %% Get the atmosphere (optional)
# The atmosphere is the full atmospheric field of ARTS as a dictionary,
# which is likely more than you wish to change.  You may change only part
# of the atmosphere by simply creating a dictionary that only contains the
# fields that you want to change.
atm = fop.get_atmosphere()

# %% Get the profile flux for the given `atm`
# Passing `atm` is optional, if not passed the operator will use the current atmosphere,
# which is the atmosphere that was set with the last call to `__call__`, or the constructor
# default if no call to `__call__` has been made.
kays = np.linspace(400, 2500, 10001)
flux, alts = fop(pyarts.arts.convert.kaycm2freq(kays), atm)

# %% Plots
fig, ax = plt.subplots(figsize=(6, 4))
ax.set_rasterization_zorder(0)
cx = ax.contourf(
    kays, alts / 1e3, flux.up.T, 50, locator=ticker.LogLocator(), zorder=-1
)
ax.set_xlabel("Kaysers [cm$^{-1}$]")
ax.set_ylabel("Altitude [km]")
ax.set_title("Upgoing spectral irradiance")
fig.colorbar(cx, label="Flux [W / m$^2$ Hz]")

fig, ax = plt.subplots(figsize=(6, 4))
ax.set_rasterization_zorder(0)
cx = ax.contourf(
    kays, alts / 1e3, flux.down.T, 50, locator=ticker.LogLocator(), zorder=-1
)
ax.set_xlabel("Kaysers [cm$^{-1}$]")
ax.set_ylabel("Altitude [km]")
ax.set_title("Downgoing spectral irradiance")
fig.colorbar(cx, label="Flux [W / m$^2$ Hz]")
