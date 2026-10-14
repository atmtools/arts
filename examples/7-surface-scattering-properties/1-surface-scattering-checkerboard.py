"""
What a checkerboard surface looks like from a high orbit.

This example drives the surface scattering framework as the lower boundary of a
complete simulation.  A map of scattering models is keyed by surface property
tags, the surface field carries one mask per tag, and the
spectral_rad_surface_agenda resolves the model of every hit point from the mask
that is set there.

The scene is an atmosphere-free Earth, so the only light reaching the sensor is
sunlight reflected by the surface and the thermal emission of that surface.
The surface is split into 45 deg x 45 deg cells and painted in a 4 x 8
checkerboard of two Lambertian types with opposite spectra:

* red   -- reflectivity 1.0 at 700 nm, falling to 0.0 at 400 nm,
* blue  -- reflectivity 0.0 at 700 nm, rising to 1.0 at 400 nm.

Both spectra live on their own frequency grid and the model interpolates them
onto f_grid (BRDF = r, emissivity = 1 - r), which is why the 650/550/400 nm
bands of the simulation land between the stored 700/600/500/400 nm points.  A
blackbody sun with its subsolar point at (11.25 N, 11.25 W) lights a 280 K
surface, and a Dirac sensor at 36 000 km altitude looks straight down at it.
Without an atmosphere there is no sky to light up, so the surface agenda uses
the SurfaceScatteringModelDirectOnly option: it reflects the sun as a
collimated beam taken straight from the suns and skips the hemispherical
quadrature that the full option integrates over.  The only thing left out is
the reflection of the cosmic background, which is nothing at these frequencies.

Two things are measured and plotted.  The spectra above one blue and one red
cell show the two reflectivity spectra imprinted on the upwelling radiance, and
how a band with r = 0 collapses to the thermal floor.  The scan over the centre
of every cell shows the checkerboard modulated by the day-night pattern of the
sun, with the night side dark at every band.

The example ends by pinning the spectra and the scan to reference arrays and by
checking that exactly the sunlit cells are bright, so it doubles as a
regression test of the framework.
"""

import os

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pyarts3 as pa

arts = pa.arts

# %% Constants

SPEED_OF_LIGHT = 299792458.0

# Ascending in frequency means descending in wavelength
BAND_WAVELENGTHS_NM = [650.0, 550.0, 400.0]
REFLECTIVITY_WAVELENGTHS_NM = [700.0, 600.0, 500.0, 400.0]
RED_REFLECTIVITY = [1.0, 0.9, 0.5, 0.0]
BLUE_REFLECTIVITY = [0.0, 0.25, 0.7, 1.0]

CELL_SIZE = 45.0
SURF_T = 280.0
TOA = 100e3
MAX_STEPSIZE = 12000.0

SUN_LAT, SUN_LON = 11.25, -11.25
SENSOR_ALT = 36e6

# Cell centres of the 4 x 8 checkerboard
CELL_LATS = [67.5, 22.5, -22.5, -67.5]
CELL_LONS = [-157.5, -112.5, -67.5, -22.5, 22.5, 67.5, 112.5, 157.5]

# Cell edges of the checkerboard, for the pcolormesh plots
LAT_EDGES = np.array([90.0, 45.0, 0.0, -45.0, -90.0])
LON_EDGES = np.arange(len(CELL_LONS) + 1) * CELL_SIZE - 180.0

# The two probe positions of the red vs blue comparison, at cell centres
BLUE_PROBE = (22.5, 22.5)
RED_PROBE = (22.5, 67.5)

# Radiance floors of the plots: values at or below are not resolved
RADIANCE_FLOOR = 1e-10
SPECTRA_FLOOR = 1e-9

RED = arts.SurfacePropertyTag("red")
BLUE = arts.SurfacePropertyTag("blue")


def nm2hz(wavelength_nm):
    """Wavelength [nm] to frequency [Hz]"""
    return SPEED_OF_LIGHT / (np.asarray(wavelength_nm, dtype=float) * 1e-9)


# %% The two surface scattering models

def reflectivity_spectrum(values):
    """Reflectivity spectrum on its own ascending frequency grid"""
    return arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[nm2hz(REFLECTIVITY_WAVELENGTHS_NM)],
        data=np.asarray(values, dtype=float),
    )


def surface_models():
    """The map of surface scattering models, keyed by the checkerboard tags"""
    models = arts.MapOfSurfaceScatteringModel()
    models.add("red", arts.LambertianSurfaceScatterer(reflectivity_spectrum(RED_REFLECTIVITY)))
    models.add("blue", arts.LambertianSurfaceScatterer(reflectivity_spectrum(BLUE_REFLECTIVITY)))
    return models


# %% The checkerboard

def red_mask(lat, lon):
    """1.0 on red cells, 0.0 on blue cells.

    Cells are [lat, lat + 45) x [lon, lon + 45) with edges at multiples of 45
    deg, and the phase puts (22.5 N, 22.5 E) on a blue cell.  The arguments are
    rounded to 1e-9 deg first, because the ray-traced hit point carries a
    1e-14 deg jitter that would otherwise flip a cell edge at random.
    """
    i = int(np.floor(round(lat, 9) / CELL_SIZE))
    j = int(np.floor(round(lon + 180.0, 9) / CELL_SIZE))
    return 1.0 if (i + j) % 2 else 0.0


# %% The workspace

def build_workspace():
    """Set up the whole simulation with the surface scattering agenda"""
    ws = pa.Workspace()

    ws.freq_grid = nm2hz(BAND_WAVELENGTHS_NM)

    # An empty atmosphere is a non-absorbing, non-scattering one
    ws.abs_species = []
    ws.abs_bands = {}
    ws.spectral_propmat_agendaAuto()
    ws.atm_fieldInit(toa=TOA)
    ws.max_stepsize = MAX_STEPSIZE

    ws.surf_fieldEarth()
    ws.surf_field["t"] = SURF_T
    ws.surf_field[RED] = red_mask
    ws.surf_field[BLUE] = lambda lat, lon: 1.0 - red_mask(lat, lon)

    ws.surface_models = surface_models()

    ws.sunBlackbody(lat=SUN_LAT, lon=SUN_LON)
    ws.suns = [ws.sun]

    ws.jac_targets = arts.JacobianTargets()

    ws.ray_path_observer_agendaSetGeometric()
    ws.spectral_rad_space_agendaSet(option="SunOrCosmicBackground")
    ws.spectral_rad_incoming_agendaSet(option="Emission")
    ws.spectral_rad_closed_surface_agendaSet(option="Blackbody")
    # The sun is the only light source worth reflecting, so the agenda takes it
    # as a beam directly from the suns instead of integrating the sky over a
    # quadrature hemisphere.  What that drops is the diffuse sky, which in an
    # empty atmosphere is just the cosmic background.
    ws.spectral_rad_surface_agendaSet(option="SurfaceScatteringModelDirectOnly")
    ws.spectral_rad_observer_agendaSet(option="Emission")

    return ws


def measure_nadir(ws, lat, lon):
    """Point a Dirac sensor at (lat, lon) straight down and simulate the spectrum"""
    ws.measurement_sensorSimple(pos=[SENSOR_ALT, lat, lon], los=[180.0, 0.0])
    ws.measurement_vecFromSensor()
    return np.array(ws.measurement_vec)


# %% The red vs blue comparison: nadir over one blue and one red cell

ws = build_workspace()
blue_spectrum = measure_nadir(ws, *BLUE_PROBE)
red_spectrum = measure_nadir(ws, *RED_PROBE)

PROBES = [
    ("blue cell (22.5 N, 22.5 E)", BLUE_PROBE, "tab:blue", blue_spectrum),
    ("red cell (22.5 N, 67.5 E)", RED_PROBE, "tab:red", red_spectrum),
]

# %% The checkerboard scan: nadir over the centre of every cell

scan = np.zeros((len(CELL_LATS), len(CELL_LONS), len(BAND_WAVELENGTHS_NM)))
for i, lat in enumerate(CELL_LATS):
    for j, lon in enumerate(CELL_LONS):
        scan[i, j, :] = measure_nadir(ws, lat, lon)


# %% Plots

def mark_probes(ax):
    """Mark the two probe positions of the red vs blue comparison on a map"""
    for label, (lat, lon), color, _ in PROBES:
        ax.plot(lon, lat, marker="*", ms=15, mfc=color, mec="white", mew=1.2, ls="none")
        ax.annotate(label.split()[0], (lon, lat), textcoords="offset points", xytext=(0, 11),
                    ha="center", fontsize=9, fontweight="bold", color="white",
                    bbox=dict(boxstyle="round,pad=0.2", fc="black", alpha=0.5))


fig_masks, ax_masks = plt.subplots(1, 2, figsize=(10, 3.6), sharex=True, sharey=True,
                                   constrained_layout=True)
red_grid = np.array([[red_mask(lat, lon) for lon in CELL_LONS] for lat in CELL_LATS])
for ax, (key, grid, color) in zip(ax_masks, (("red", red_grid, "#b2182b"),
                                             ("blue", 1.0 - red_grid, "#2166ac"))):
    ax.pcolormesh(LON_EDGES, LAT_EDGES, grid, cmap=ListedColormap(["white", color]),
                  vmin=0.0, vmax=1.0)
    ax.set_title(f"{key} cells")
    ax.set_xlabel("Longitude [deg]")
    ax.set_xticks(np.arange(-180.0, 180.0, 45.0))
    ax.set_yticks(np.arange(-90.0, 90.0, 45.0))
    ax.set_box_aspect(0.5)
    mark_probes(ax)
ax_masks[0].set_ylabel("Latitude [deg]")
fig_masks.suptitle("The checkerboard surface masks")

fig_refl, ax_refl = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
for label, values, color in (("red", RED_REFLECTIVITY, "tab:red"),
                             ("blue", BLUE_REFLECTIVITY, "tab:blue")):
    ax_refl.plot(REFLECTIVITY_WAVELENGTHS_NM[::-1], values[::-1], "o-",
                 ms=7, lw=1.5, color=color, label=f"{label} spectrum")
    # Linear interpolation in frequency, exactly what the model does onto f_grid
    r_bands = np.interp(nm2hz(BAND_WAVELENGTHS_NM), nm2hz(REFLECTIVITY_WAVELENGTHS_NM), values)
    ax_refl.plot(BAND_WAVELENGTHS_NM[::-1], r_bands[::-1], "s", ms=9, mfc="white", mec=color,
                 mew=1.5, ls="none", label=f"{label} on f_grid")
    for w, r in zip(BAND_WAVELENGTHS_NM[::-1], r_bands[::-1]):
        ax_refl.annotate(f"{r:.2f}", (w, r), textcoords="offset points", xytext=(0, 9),
                         ha="center", fontsize=7, color=color)
for wl in BAND_WAVELENGTHS_NM:
    ax_refl.axvline(wl, color="0.6", ls=":", lw=0.8)
ax_refl.set_xlabel("Wavelength [nm]")
ax_refl.set_ylabel("Reflectivity r(λ) [-]")
ax_refl.set_xlim(380, 720)
ax_refl.set_ylim(0, 1.12)
ax_refl.set_title("Surface reflectivity spectra, sampled and interpolated onto f_grid")
ax_refl.grid(True, alpha=0.3)
ax_refl.legend(ncols=2, framealpha=0.9, edgecolor="0.6", fontsize=8)

# Floor the radiance at RADIANCE_FLOOR so that the invisible thermal emission
# does not stretch the log colour scale decades beyond the reflected sunlight
log_scan = np.log10(np.maximum(scan, RADIANCE_FLOOR))
fig_scan, ax_scan = plt.subplots(1, 3, figsize=(12.5, 4.0), sharex=True, sharey=True,
                                 constrained_layout=True)
for ax, k in zip(ax_scan, (2, 1, 0)):  # panels in ascending wavelength
    mesh = ax.pcolormesh(LON_EDGES, LAT_EDGES, log_scan[:, :, k], cmap="magma",
                         vmin=np.log10(RADIANCE_FLOOR), vmax=log_scan.max())
    ax.plot(SUN_LON, SUN_LAT, marker="o", ms=10, mfc="#ff7f0e", mec="white", mew=1.2, ls="none")
    ax.set_title(f"{BAND_WAVELENGTHS_NM[k]:.0f} nm radiance")
    ax.set_xlabel("Longitude [deg]")
    ax.set_xticks(np.arange(-180.0, 180.0, 45.0))
    ax.set_yticks(np.arange(-90.0, 90.0, 45.0))
    ax.set_box_aspect(0.5)
    mark_probes(ax)
ax_scan[0].annotate("subsolar", (SUN_LON, SUN_LAT), textcoords="offset points", xytext=(12, -16),
                    fontsize=9, fontweight="bold", color="gold",
                    bbox=dict(boxstyle="round,pad=0.2", fc="black", alpha=0.5))
ax_scan[0].set_ylabel("Latitude [deg]")
fig_scan.colorbar(mesh, ax=ax_scan.tolist(), pad=0.015,
                  label=r"log$_{10}$ radiance [W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$]")
fig_scan.suptitle("Nadir radiance over the checkerboard at 36 000 km")

fig_spectra, ax_spectra = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
wl = BAND_WAVELENGTHS_NM[::-1]  # ascending wavelength
for label, _, color, spectrum in PROBES:
    vals = spectrum[::-1]
    ax_spectra.semilogy(wl, np.where(vals > SPECTRA_FLOOR, vals, np.nan), "o-",
                        ms=7, lw=1.5, label=label, color=color)
    for w, v in zip(wl, vals):
        if v > SPECTRA_FLOOR:
            ax_spectra.annotate(f"{v:.2e}", (w, v), textcoords="offset points", xytext=(0, 10),
                                ha="center", fontsize=7, color=color)
ax_spectra.annotate("r = 0 → dark", (400.0, SPECTRA_FLOOR * 1.5), ha="center",
                    fontsize=7, color="0.45")
ax_spectra.set_xlabel("Wavelength [nm]")
ax_spectra.set_ylabel("Radiance [W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$]")
ax_spectra.set_ylim(SPECTRA_FLOOR, 3e-8)
ax_spectra.set_title("Red vs blue surface, nadir at 36 000 km")
ax_spectra.grid(True, which="both", alpha=0.3)
ax_spectra.legend(framealpha=0.9, edgecolor="0.6")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% References

ref_blue = [3.0520217387884524e-09, 9.773891841217284e-09, 1.0126225482288864e-08]
ref_red = [7.659689944605077e-09, 4.688325511370435e-09, 1.005756739358586e-61]
ref_scan = [
    [
        [3.1045589331336524e-42, 1.7911284030340405e-47, 1.005756739358586e-61],
        [3.9428034091852686e-10, 1.2626559497527454e-09, 1.3081727382976113e-09],
        [1.1844015213713968e-08, 7.2494577567888405e-09, 1.005756739358586e-61],
        [2.0209971276512035e-09, 6.472105714723686e-09, 6.7054150871737055e-09],
        [1.4999280630627336e-08, 9.180725400289066e-09, 1.005756739358586e-61],
        [9.340529343350676e-10, 2.9912409332269738e-09, 3.099070529302691e-09],
        [3.1045589331336524e-42, 1.7911284030340405e-47, 1.005756739358586e-61],
        [5.950404621839493e-41, 3.4667001349045923e-47, 0.0],
    ],
    [
        [5.950404621839493e-41, 3.4667001349045923e-47, 0.0],
        [3.1045589331336524e-42, 1.7911284030340405e-47, 1.005756739358586e-61],
        [2.130531202964453e-09, 6.822881133991974e-09, 7.068835416236065e-09],
        [2.9352759331153067e-08, 1.7966169831493977e-08, 1.005756739358586e-61],
        [3.0520217387884524e-09, 9.773891841217284e-09, 1.0126225482288864e-08],
        [7.659689944605077e-09, 4.688325511370435e-09, 1.005756739358586e-61],
        [5.950404621839493e-41, 3.4667001349045923e-47, 0.0],
        [3.1045589331336524e-42, 1.7911284030340405e-47, 1.005756739358586e-61],
    ],
    [
        [3.1045589331336524e-42, 1.7911284030340405e-47, 1.005756739358586e-61],
        [5.950404621839493e-41, 3.4667001349045923e-47, 0.0],
        [1.3062727697291503e-08, 7.99540450778057e-09, 1.005756739358586e-61],
        [3.000356403250418e-09, 9.608437121458504e-09, 9.95480637651148e-09],
        [2.0680285802667647e-08, 1.2657942059308443e-08, 1.005756739358586e-61],
        [3.762281953541073e-10, 1.2048451825468122e-09, 1.2482779826013501e-09],
        [3.1045589331336524e-42, 1.7911284030340405e-47, 1.005756739358586e-61],
        [5.950404621839493e-41, 3.4667001349045923e-47, 0.0],
    ],
    [
        [5.950404621839493e-41, 3.4667001349045923e-47, 0.0],
        [3.1045589331336524e-42, 1.7911284030340405e-47, 1.005756739358586e-61],
        [1.0408450214738082e-10, 3.3332353220370097e-10, 3.4533932853786043e-10],
        [5.723176787246221e-09, 3.5030289648510007e-09, 1.005756739358586e-61],
        [4.857581029128194e-10, 1.5556072548649289e-09, 1.6116844835766866e-09],
        [3.1045589331336524e-42, 1.7911284030340405e-47, 1.005756739358586e-61],
        [5.950404621839493e-41, 3.4667001349045923e-47, 0.0],
        [3.1045589331336524e-42, 1.7911284030340405e-47, 1.005756739358586e-61],
    ],
]

# The sub-floor thermal terms are only loosely pinned (atol dwarfs them); the
# reflected sunlight, which is what the plots show, is pinned to relative
# machine precision
assert np.allclose(blue_spectrum, ref_blue, rtol=1e-8, atol=1e-30), \
    f"Mismatch to blue cell reference.\nReference:   {ref_blue}\nSimulation:  {blue_spectrum}"

assert np.allclose(red_spectrum, ref_red, rtol=1e-8, atol=1e-30), \
    f"Mismatch to red cell reference.\nReference:   {ref_red}\nSimulation:  {red_spectrum}"

assert np.allclose(scan, ref_scan, rtol=1e-8, atol=1e-30), \
    "Mismatch to the scan reference"

# The scan is the checkerboard modulated by the sun: every sunlit cell glows
# with reflected sunlight above the floor, every night-side cell stays below
# it.  The 1e-12 tolerance keeps the two cells sitting exactly on the horizon
# (cos of the incidence angle is zero to machine precision) unambiguously dark
for i, lat in enumerate(CELL_LATS):
    for j, lon in enumerate(CELL_LONS):
        cos_inc = (np.sin(np.deg2rad(SUN_LAT)) * np.sin(np.deg2rad(lat)) +
                   np.cos(np.deg2rad(SUN_LAT)) * np.cos(np.deg2rad(lat)) *
                   np.cos(np.deg2rad(lon - SUN_LON)))
        sunlit = cos_inc > 1e-12
        bright = bool(np.any(scan[i, j] > RADIANCE_FLOOR))
        assert bright == sunlit, f"Cell ({lat}, {lon}) sunlit={sunlit}, scan={scan[i, j]}"
