# %%
"""
Simulation of a photo of Earth at 500 GHz from the Moon using the ARTS camera model.

This script uses the thin-lens camera sensor model to simulate the radiance
observed by a camera at the Moon looking at Earth, creating a 2D image.
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pa

# 1. Simulation Parameters
freq_hz = 500e9  # 500 GHz
moon_dist = 384400e3  # Average Moon-Earth distance in meters
earth_radius = 6371e3  # Earth radius in meters

# Angular radius of Earth as seen from Moon distance
earth_angular_radius_deg = np.rad2deg(np.arctan2(earth_radius, moon_dist))  # ~0.95 degrees

# Camera parameters
n_h = 40  # pixel rows
n_w = 40  # pixel columns
focal_length = 0.050  # 50 mm lens
aperture_diameter = 0.025  # 25 mm aperture (f/2)
focus_distance = moon_dist  # Focus at Earth distance

# Choose CCD size so the field of view covers ~1.5x the Earth angular diameter.
# For focus at near-infinity: image_distance ≈ focal_length
# => ccd_half ≈ focal_length * tan(fov_half)
fov_half_rad = np.deg2rad(earth_angular_radius_deg * 1.5)
ccd_half = focal_length * np.tan(fov_half_rad)
ccd_h = 2.0 * ccd_half  # symmetric CCD
ccd_w = 2.0 * ccd_half

print(f"Earth angular radius from Moon: {earth_angular_radius_deg:.4f} degrees")
print(f"Camera FOV half-angle: {np.rad2deg(fov_half_rad):.4f} degrees")
print(f"CCD size: {ccd_h*1e3:.2f} x {ccd_w*1e3:.2f} mm")
print(f"Total pixels: {n_h} x {n_w} = {n_h * n_w}")

# 2. ARTS Workspace Setup
ws = pa.workspace.Workspace()

# Set frequency grid (single frequency for this simulation)
ws.freq_grid = np.array([freq_hz])

# Select absorption species
ws.abs_speciesSet(
    species=["O2-PWR98", "H2O-PWR98"]
)

# Download and read spectral line data
ws.ReadCatalogData()

# Apply frequency cutoff to speed up calculations
cutoff = pa.arts.convert.kaycm2freq(25)
for band in ws.abs_bands:
    ws.abs_bands[band].cutoff = "ByLine"
    ws.abs_bands[band].cutoff_value = cutoff

# Automatically set up the methods to compute absorption coefficients
ws.spectral_propmat_agendaAuto()

# Set up a simple atmosphere
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pa.arts.SurfaceKey("t")] = 288.0

# Load a standard atmosphere profile
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)

# 3. Build camera sensor using the new model
# Observer at Moon distance, looking down at Earth
pos = [moon_dist, 0.0, 0.0]    # altitude = Moon distance
los = [180.0, 0.0]              # zenith=180 = looking straight down

ws.measurement_sensorsInit()
ws.measurement_sensorsAddCamera(
    pos=pos,
    los=los,
    n_h=n_h,
    n_w=n_w,
    ccd_h=ccd_h,
    ccd_w=ccd_w,
    focal_length=focal_length,
    aperture_diameter=aperture_diameter,
    focus_distance=focus_distance,
)

print(f"Sensor elements: {len(ws.measurement_sensors[0].elements)}")

# 4. Run forward model
ws.spectral_rad_transform_operatorSet(option="Tb")
ws.max_stepsize = 1000.0
ws.ray_path_observer_agendaSetGeometric()

n_total = len(ws.measurement_sensors[0].elements)
print(f"Running forward model for {n_total} sensor elements ({n_h}x{n_w} pixels x {len(ws.freq_grid)} freq)...")

ws.measurement_vecFromSensor(kernel="High Performance")

# 5. Give measurement_vec the structure of the camera's meta and extract the image
camera_gf = ws.measurement_sensors[0].reshape(ws.measurement_vec)

# camera_gf.data is shaped (n_h, n_w, n_freq); grids give angular offsets
image_data = np.array(camera_gf.data)[:, :, 0]
row_offsets = np.array(camera_gf.grids[0])  # row angular offsets (degrees)
col_offsets = np.array(camera_gf.grids[1])  # col angular offsets (degrees)

print(f"Camera image shape: {image_data.shape}")
print(f"Row offset range: [{row_offsets[0]:.4f}, {row_offsets[-1]:.4f}] deg")
print(f"Col offset range: [{col_offsets[0]:.4f}, {col_offsets[-1]:.4f}] deg")

# 6. Visualization
col_grid, row_grid = np.meshgrid(col_offsets, row_offsets)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Subplot 1: Scatter plot using the angular offset grids from meta
sc = ax1.scatter(
    col_grid.ravel(),
    row_grid.ravel(),
    c=image_data.ravel(),
    s=80,
    cmap="inferno",
    edgecolors="none",
)
# Draw Earth limb circle
theta_circ = np.linspace(0, 2 * np.pi, 200)
ax1.plot(
    earth_angular_radius_deg * np.cos(theta_circ),
    earth_angular_radius_deg * np.sin(theta_circ),
    "c--",
    lw=1.5,
    label=f"Earth limb ({earth_angular_radius_deg:.2f}°)",
)
ax1.set_xlabel("Column angular offset (degrees)")
ax1.set_ylabel("Row angular offset (degrees)")
ax1.set_title(f"Earth from Moon at {freq_hz/1e9:.0f} GHz — Camera image")
ax1.set_aspect("equal", "box")
ax1.legend(loc="upper right")
plt.colorbar(sc, ax=ax1, label="Brightness Temperature (K)")
ax1.grid(True, linestyle="--", alpha=0.3)

# Subplot 2: imshow using the angular offset grids as extent
extent = [col_offsets[0], col_offsets[-1], row_offsets[0], row_offsets[-1]]
im2 = ax2.imshow(
    image_data,
    origin="lower",
    cmap="inferno",
    aspect="equal",
    extent=extent,
)
ax2.set_xlabel("Column angular offset (degrees)")
ax2.set_ylabel("Row angular offset (degrees)")
ax2.set_title(f"Earth from Moon at {freq_hz/1e9:.0f} GHz — Pixel grid")
plt.colorbar(im2, ax=ax2, label="Brightness Temperature (K)")

plt.tight_layout()

if "ARTS_HEADLESS" not in os.environ:
    plt.show()
