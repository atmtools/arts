# %%
"""
Calculate the clear sky upwelling longwave radiance at the top of the atmosphere
and, for each frequency, the height at which the optical depth along
the downlooking ray reaches unity (the "emission height" or tau=1 level).

The atmosphere is a standard tropical atmosphere, containing water vapor
(with continuum), carbon dioxide, and ozone as trace gases.
"""

from os import environ
import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pa
from scipy.constants import speed_of_light

# Number of frequencies to use
nf = 1000

# Minimum total optical depth required to define an emission height
TAU_TOTAL_MIN = 0.0

# Radiance from ARTS is per Hz; convert to per Kayser (cm^-1) via d(nu) = c d(kayser)
c_cm_s = speed_of_light * 100.0  # speed of light in cm/s

# Download ARTS catalogs if they are not already present
pa.data.download()

# ARTS workspace
ws = pa.workspace.Workspace()

# Set up frequency grid
kayser_grid = np.linspace(1, 2000, nf)  # in Kayser (cm^-1)
ws.freq_grid = pa.arts.convert.kaycm2freq(kayser_grid)  # in Hz

# Select absorption species and continuum model
# This example uses a reduced set of species to speed up the calculation.
# Use the second line for a more realistic setup.
ws.abs_speciesSet(
    species=["H2O-161", "H2O-ForeignContCKDMT400", "H2O-SelfContCKDMT400", "CO2-626"]
)
# ws.abs_speciesSet(
#     species=["H2O", "H2O-ForeignContCKDMT400", "H2O-SelfContCKDMT400", "CO2", "O3"]
# )

# Read spectral line data from ARTS catalog
ws.ReadCatalogData()

# Apply a frequency cutoff. To be consistent with the CKD water vapor continuum,
# a cutoff of 25 Kayser is necessary. We set it here for all species, because it
# also speeds up the calculation.
cutoff = pa.arts.convert.kaycm2freq(25)
for band in ws.abs_bands:
    ws.abs_bands[band].cutoff = "ByLine"
    ws.abs_bands[band].cutoff_value = cutoff

# Remove 90% of the lines to speed up the calculation
ws.abs_bands.keep_hitran_s(approximate_percentile=90)

# Automatically set up the methods to compute absorption coefficients
ws.spectral_propmat_agendaAuto()

# Set up a simple atmosphere
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pa.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1)

# Set up geometry of observation: downlooking from the top of the atmosphere
pos = [100e3, 0, 0]

# Single angle from the Gauss-Legendre quadrature
los = [
    180.0 - np.degrees(np.arccos(1 / np.sqrt(3))),  # ~125°
    0.0,
]
ws.ray_pathGeometric(pos=pos, los=los, max_stepsize=1000.0)

# %% Compute the tau=1 height per frequency

# Build the per-level atmospheric state and propagation matrices along the path.
# These are the building blocks used internally by spectral_radClearskyEmission;
# we stop before the emission solve and keep the per-level absorption data.
ws.ray_pointBackground()
ws.atm_pathFromPath()
ws.freq_grid_pathFromPath()
ws.spectral_propmat_pathFromPath()

# Absorption coefficient [1/m] at every path point and frequency.
K = np.array(
    [
        [ws.spectral_propmat_path[ip][f][0] for f in range(nf)]
        for ip in range(len(ws.ray_path))
    ]
)  # shape [np, nf]

# Altitude of every path point [m]
alt = np.array([p.pos[0] for p in ws.ray_path])  # shape [np]

# Per-layer path length [m] (length np-1)
r = ws.ray_path.distances(ws.surf_field.ellipsoid)

# Per-layer optical depth via the trapezoidal rule, matching the layer
# transmittance used internally by spectral_tramat_pathFromPath.
dtau = 0.5 * (K[:-1] + K[1:]) * r[:, None]  # shape [np-1, nf]

# Cumulative optical depth from the observer (TOA) downwards.
tau_cum = np.zeros((len(alt), nf))
tau_cum[1:] = np.cumsum(dtau, axis=0)  # tau_cum[0] == 0 (observer level)

# Per-frequency total optical depth from TOA down to the surface.
tau_total = tau_cum[-1, :]  # shape [nf]

# For each frequency, find the first level where optical depth
# crosses 1 and linearly interpolate the altitude at the crossing.
emission_height = np.full(nf, np.nan)
for f in range(nf):
    if tau_total[f] < TAU_TOTAL_MIN:
        continue
    t = tau_cum[:, f]
    # ray_path is ordered from TOA (index 0) to surface (index -1), so
    # tau increases with index; searchsorted gives the first
    # index where t >= 1.
    idx = np.searchsorted(t, 1.0)
    idx = np.clip(idx, 0, len(alt) - 1)
    x = np.interp(1.0, [t[idx - 1], t[idx]], [alt[idx - 1], alt[idx]])
    emission_height[f] = x

# %% Calculate the clear sky upwelling longwave radiance at the top of the atmosphere
ws.spectral_radClearskyEmission()

# %% Show results


# Spectra are noisy at the line-by-line resolution, so plot the full
# resolution data in light gray and overlay a smoothed curve on top.
def smooth(y, window=25):
    if window <= 1:
        return y
    kernel = np.ones(window) / window
    ypad = np.pad(y, window // 2, mode="edge")
    return np.convolve(ypad, kernel, mode="valid")


# 25 Kayser smoothing window
window_size = int(round(25 / np.diff(kayser_grid[:2]).item()))
spectral_rad_i = np.array([ws.spectral_rad[f][0] for f in range(nf)])
emission_height_smooth = smooth(emission_height / 1e3, window=window_size)[:nf]
spectral_rad_smooth = smooth(spectral_rad_i, window=window_size)[:nf]
olr_single = np.trapezoid(spectral_rad_i * c_cm_s, kayser_grid) * np.pi

fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8, 8))

ax_top.plot(kayser_grid, spectral_rad_i * c_cm_s, color="lightgray", linewidth=1.0)
ax_top.plot(kayser_grid, spectral_rad_smooth * c_cm_s, color="black", linewidth=1.5)
ax_top.set_xlabel("Wavenumber / cm$^{-1}$")
ax_top.set_ylabel("Spectral radiance / W m$^{-2}$ sr$^{-1}$ cm$^{-1}$")
ax_top.set_title(
    f"Clear sky outgoing radiance los={los[0]:.1f}° OLR={olr_single:.2f} Wm$^{{-2}}$"
)

ax_bot.plot(kayser_grid, emission_height / 1e3, color="lightgray", linewidth=1.0)
ax_bot.plot(kayser_grid, emission_height_smooth, color="black", linewidth=1.5)
ax_bot.set_xlabel("Wavenumber / cm$^{-1}$")
ax_bot.set_ylabel("Emission height / km")
ax_bot.set_title(r"Height where optical depth reaches $\tau = 1$")
ax_bot.set_ylim(0, 30)
fig.tight_layout()

if "ARTS_HEADLESS" not in environ:
    plt.show()
