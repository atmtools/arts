import pyarts
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

# Initialize the operator
toa = pyarts.recipe.SpectralAtmosphericFlux()

# Get the atmosphere (optional)
# The atmosphere is the full atmospheric field of ARTS as a dictionary,
# which is likely more than you wish to change.  You may change only part
# of the atmosphere by simply creating a dictionary that only contains the
# fields that you want to change.
atm = toa.get_atmosphere()

# Get the top-of-the-atmosphere flux for the given atm
# Passing "atm" is optional, if not passed the operator will use the current atmosphere,
# which is the atmosphere that was set with the last call to __call__, or the constructor
# default if no call to __call__ has been made.
kays = np.linspace(400, 2500, 10001)
flux, alts = toa(pyarts.arts.convert.kaycm2freq(kays), atm)

# Plot the results
plt.figure(1, figsize=(6, 4))
plt.clf()
plt.contourf(kays, alts/1e3, flux.up.T, 50, locator=ticker.LogLocator())
plt.xlabel("Kaysers [cm-1")
plt.ylabel("Altitude [km")
cm = plt.colorbar()
cm.set_label("Flux [W / m$^2$ Hz]")
plt.title("Upgoing spectral irradiance")

plt.figure(1, figsize=(6, 4))
plt.clf()
plt.contourf(kays, alts/1e3, flux.down.T, 50, locator=ticker.LogLocator())
plt.xlabel("Kaysers [cm-1")
plt.ylabel("Altitude [km")
cm = plt.colorbar()
cm.set_label("Flux [W / m$^2$ Hz]")
plt.title("Downgoing spectral irradiance")
