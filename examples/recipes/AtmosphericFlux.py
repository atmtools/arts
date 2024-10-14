import pyarts
import matplotlib.pyplot as plt

# Initialize the operator
toa = pyarts.recipe.AtmosphericFlux()

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
solar, thermal, altitude = toa(atmospheric_profile=atm)

# Plot the results
plt.figure(1, figsize=(8, 6))
plt.clf()
plt.subplot(121)
plt.plot(solar.up, altitude/1e3)
plt.plot(solar.down, altitude/1e3)
plt.legend(["up", "down"])
plt.ylabel("Altitude [km]")
plt.xlabel("Flux [W / m$^2$]")
plt.title("Solar flux")

plt.subplot(122)
plt.plot(thermal.up, altitude/1e3)
plt.plot(thermal.down, altitude/1e3)
plt.legend(["up", "down"])
plt.ylabel("Altitude [km]")
plt.xlabel("Flux [W / m$^2$]")
plt.title("Thermal flux")
plt.show()
