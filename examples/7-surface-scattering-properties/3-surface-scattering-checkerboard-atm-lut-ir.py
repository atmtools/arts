"""
What a checkerboard surface looks like from a low orbit through an atmosphere, in the thermal infrared.

This is the infrared counterpart of ``2-surface-scattering-checkerboard-atm-lut.py``:
the same 4 x 8 checkerboard of red and blue Lambertian cells, the same Dirac
sensor looking nadir, and the same tropical clearsky atmosphere up to 100 km
with the absorption taken from a lookup table -- but the sun and the Rayleigh
scattering are gone, and the measurement bands move from three visible
wavelengths to two window wavenumbers, 895 and 1145 /cm.

What the sensor sees is now thermal emission rather than reflected sunlight.
The surface emits its own Kirchhoff radiance, emissivity = 1 - r times the
blackbody at the local surface temperature, upwards into the atmosphere, where
it is absorbed on the way up while the air itself glows into the path; the
diffuse sky coming back down is still reflected at the surface over the same
coarse quadrature hemisphere.  The LambertianSurfaceScatterer carries the
emissivity implicitly: its reflectivity r(f) fixes emissivity(f) = 1 - r(f),
and the SurfaceScatteringModelDiffuseOnly surface agenda adds each channel's
thermal emission term together with the reflected sky.  With no sun in the
game there is no direct beam to count separately, and the incoming and
observer agendas are the plain Emission options -- clearsky emission along the
path, no scattering source -- so no hand-written agendas are needed.

The emissivity spectra mirror the reflectivity trend of the solar example: the
red cells emit strongly at 895 /cm (emissivity 0.99) and weakly at 1145 /cm
(0.90), the blue cells exactly the other way (0.90 / 0.99).  The absorbing
species are H2O, O3 and CO2, and the line window covers both bands, so the
9.6 micron ozone band at about 1043 /cm sits right between the measurements.
That is what separates the two bands: 895 /cm is a clean window that sees all
the way down to the surface, while 1145 /cm couples to the ozone layer and
sees further up the colder atmosphere -- the cell-to-cell emissivity contrast
is strong at the first band and damped at the second.

The surface temperature is no longer uniform either.  It falls quadratically
from 300 K at the equator to 250 K at the poles, and the red cells carry a
+10 K offset, so the scan over the cell centres shows the checkerboard riding
on top of a latitudinal brightness gradient instead of the day-night pattern
of the solar example.

Everything is plotted on linear scales -- in the thermal IR the radiances
live in one narrow band, so no log scale is needed to see the structure --
and each plot comes in two units: the spectral radiance itself, and the
Planck brightness temperature at the band frequency, the natural radiance
unit for a sounder.  The brightness temperature makes the physics directly
readable: the surface cells glow near their own 250-310 K at 895 /cm, while
the 1145 /cm panel sits systematically colder because it couples to the
ozone layer.

The example ends by pinning the spectra and the scan to reference arrays and
by checking the physical ordering -- red brighter than blue at both bands, the
red/blue contrast stronger at 895 /cm than at 1145 /cm, and the equatorial
rows brighter than the polar ones -- so it doubles as a regression test of the
framework.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pa

arts = pa.arts

# Make the arts-xml-data / arts-cat-data catalogs findable: reuses an existing
# checkout in ~/.cache/arts (or downloads one).  No-op when ARTS_DATA_PATH or
# ARTS_INCLUDE_PATH is already set in the environment, as under ctest.
pa.data.download()

# %% Constants

SPEED_OF_LIGHT = 299792458.0

# Ascending in frequency means ascending in wavenumber
BAND_WAVENUMBERS_CM = [895.0, 1145.0]
RED_EMISSIVITY = [0.99, 0.90]
BLUE_EMISSIVITY = [0.90, 0.99]

CELL_SIZE = 45.0
T_EQUATOR, T_POLE = 300.0, 250.0
RED_T_OFFSET = 10.0
TOA = 100e3
MAX_STEPSIZE = 9600.0

SENSOR_ALT = 800e3

# Cell centres of the 4 x 8 checkerboard
CELL_LATS = [67.5, 22.5, -22.5, -67.5]
CELL_LONS = [-157.5, -112.5, -67.5, -22.5, 22.5, 67.5, 112.5, 157.5]

# Cell edges of the checkerboard, for the pcolormesh plots
LAT_EDGES = np.array([90.0, 45.0, 0.0, -45.0, -90.0])
LON_EDGES = np.arange(len(CELL_LONS) + 1) * CELL_SIZE - 180.0

# The two probe positions of the red vs blue comparison, at cell centres
BLUE_PROBE = (22.5, 22.5)
RED_PROBE = (22.5, 67.5)

# Absorbing species and the line-selection window covering the two measurement
# bands (26.81 and 34.31 THz) with a margin on each side.  O3 is what makes
# the window inhomogeneous: its 9.6 micron band at ~1043 /cm falls between the
# two bands, H2O fills the window with lines, and CO2 is essentially
# transparent here (its 15 micron band at 667 /cm sits below the window).
ABS_SPECIES = ["H2O-161", "O3-666", "CO2-626"]
BAND_FMIN = 26.0e12
BAND_FMAX = 35.5e12

# The LUT grids: altitude through the whole atmosphere, and the temperature
# and humidity perturbations the line-by-line model is precomputed for.  The
# humidity axis multiplies the VMR of water_affected_species, so it only needs
# to span the factor the tropical profile can deviate by.
LUT_ALT_GRID = np.linspace(0.0, TOA, 51)
LUT_T_PERTURBATION = np.linspace(-30.0, 30.0, 5)
LUT_WATER_PERTURBATION = np.logspace(-1.0, 1.0, 5)
LUT_WATER_AFFECTED_SPECIES = ["H2O"]

# Quadrature steps of the diffuse sky integral at the surface hits
ZENITH_STEP = 45.0
AZIMUTH_STEP = 90.0

RED = arts.SurfacePropertyTag("red")
BLUE = arts.SurfacePropertyTag("blue")


def wn2hz(wavenumber_cm):
    """Wavenumber [1/cm] to frequency [Hz]"""
    return SPEED_OF_LIGHT * (np.asarray(wavenumber_cm, dtype=float) * 100.0)


def radiance_to_tb(radiance, freq):
    """Planck brightness temperature [K] of a spectral radiance
    [W m-2 sr-1 Hz-1] at a frequency [Hz] -- the inverse Planck function at
    the band frequency, the natural radiance unit in the thermal IR."""
    h = 6.62607015e-34
    k = 1.380649e-23
    return h * freq / (k * np.log1p(2 * h * freq ** 3 / (SPEED_OF_LIGHT ** 2 * radiance)))


# %% The two surface scattering models

def emissivity_spectrum(values):
    """Emissivity spectrum on its own ascending frequency grid.

    The Lambertian model is built from the reflectivity and implements
    Kirchhoff's law, emissivity = 1 - r, so the emissivity is passed as 1 - r.
    """
    return arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[wn2hz(BAND_WAVENUMBERS_CM)],
        data=1.0 - np.asarray(values, dtype=float),
    )


def surface_models():
    """The map of surface scattering models, keyed by the checkerboard tags"""
    models = arts.MapOfSurfaceScatteringModel()
    models.add("red", arts.LambertianSurfaceScatterer(emissivity_spectrum(RED_EMISSIVITY)))
    models.add("blue", arts.LambertianSurfaceScatterer(emissivity_spectrum(BLUE_EMISSIVITY)))
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


def surface_temperature(lat, lon):
    """Quadratic latitudinal temperature, 300 K equator to 250 K poles, plus
    the +10 K offset on the red cells"""
    t = T_EQUATOR - (T_EQUATOR - T_POLE) * (lat / 90.0) ** 2
    return t + RED_T_OFFSET * red_mask(lat, lon)


# %% The quadrature of the diffuse sky

def quadrature_grids():
    """Cosine-hemisphere quadrature grids and weights.

    The weights are raw multipliers in ARTS, so the ``cos * sin`` of the solid
    angle is folded into them: ``w = 2 cos(za) sin(za) dza`` for the zenith and
    ``w = daz / 2 pi`` for the azimuth, giving ``sum(w) ~ 1``.
    """
    za = np.arange(0.0, 90.0 + 1e-9, ZENITH_STEP)
    aa = np.arange(0.0, 360.0, AZIMUTH_STEP)
    za_w = 2.0 * np.cos(np.deg2rad(za)) * np.sin(np.deg2rad(za)) * np.deg2rad(ZENITH_STEP)
    aa_w = np.full(len(aa), np.deg2rad(AZIMUTH_STEP) / (2 * np.pi))
    return arts.ZenGrid(za), arts.AziGrid(aa), arts.Vector(za_w), arts.Vector(aa_w)


# %% The workspace

def build_workspace():
    """Set up the whole simulation with the atmosphere and the LUT agenda"""
    ws = pa.Workspace()

    ws.freq_grid = wn2hz(BAND_WAVENUMBERS_CM)

    # Absorbing species and the line bands covering the two measurement bands
    ws.abs_speciesSet(species=ABS_SPECIES)
    ws.ReadCatalogData()
    ws.abs_bandsSelectFrequencyByLine(fmin=BAND_FMIN, fmax=BAND_FMAX)

    # Tropical clearsky atmosphere up to TOA
    ws.atm_fieldRead(toa=TOA, basename="planets/Earth/afgl/tropical/", missing_is_zero=1)

    # The absorption lookup table, always recalculated here at the start of the
    # script so that the example carries no cached state.  This precomputes the
    # line-by-line absorption on the altitude x temperature x humidity grid of
    # the constants above; the ray tracing below only interpolates it.
    ws.abs_lookup_dataCalc(lat=0.0,
                           lon=0.0,
                           alt_grid=LUT_ALT_GRID,
                           temperature_perturbation=LUT_T_PERTURBATION,
                           water_perturbation=LUT_WATER_PERTURBATION,
                           water_affected_species=LUT_WATER_AFFECTED_SPECIES)

    # Absorption from the lookup table.  The frequency interpolation order is 0
    # because the LUT frequency grid is exactly f_grid; pressure, temperature
    # and humidity interpolate at high order between the coarse LUT nodes.
    ws.spectral_propmat_agendaAuto(use_abs_lookup_data=1,
                                   f_interp_order=0,
                                   p_interp_order=5,
                                   t_interp_order=4,
                                   water_interp_order=4)

    # The agendas require an (empty) jac_targets
    ws.jac_targets = arts.JacobianTargets()

    ws.surf_fieldEarth()
    ws.surf_field["t"] = surface_temperature
    ws.surf_field[RED] = red_mask
    ws.surf_field[BLUE] = lambda lat, lon: 1.0 - red_mask(lat, lon)

    ws.surface_models = surface_models()

    # No sun in the game, but the surface scattering methods take *suns* as a
    # hard input -- an empty list keeps the quadrature integration sun-free.
    ws.suns = arts.ArrayOfSun()

    ws.zen_grid, ws.az_grid, ws.zen_grid_weights, ws.az_grid_weights = quadrature_grids()

    ws.max_stepsize = MAX_STEPSIZE

    ws.ray_path_observer_agendaSetGeometric()
    # No sun: space is just the (negligible at these frequencies) cosmic
    # background, and both the incoming and the observer paths are plain
    # clearsky emission -- no scattering source, no sun paths to recompute.
    ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")
    ws.spectral_rad_closed_surface_agendaSet(option="Blackbody")
    ws.spectral_rad_incoming_agendaSet(option="Emission")
    # The diffuse-only option integrates the sky over the quadrature hemisphere
    # and adds each channel's thermal emission (emissivity * blackbody at the
    # local surface temperature); with no sun there is no direct part to chain.
    ws.spectral_rad_surface_agendaSet(option="SurfaceScatteringModelDiffuseOnly")
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

scan = np.zeros((len(CELL_LATS), len(CELL_LONS), len(BAND_WAVENUMBERS_CM)))
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


# The radiance maps on linear colour scales -- in the IR everything glows
# thermally, so unlike the solar example no floor is needed to tame a log scale
fig_scan, ax_scan = plt.subplots(1, 2, figsize=(9.0, 4.0), sharex=True, sharey=True,
                                 constrained_layout=True)
for ax, k in zip(ax_scan, range(len(BAND_WAVENUMBERS_CM))):
    mesh = ax.pcolormesh(LON_EDGES, LAT_EDGES, scan[:, :, k], cmap="magma",
                         vmin=scan.min(), vmax=scan.max())
    ax.set_title(f"{BAND_WAVENUMBERS_CM[k]:.0f} /cm radiance")
    ax.set_xlabel("Longitude [deg]")
    ax.set_xticks(np.arange(-180.0, 180.0, 45.0))
    ax.set_yticks(np.arange(-90.0, 90.0, 45.0))
    ax.set_box_aspect(0.5)
    mark_probes(ax)
ax_scan[0].set_ylabel("Latitude [deg]")
fig_scan.colorbar(mesh, ax=ax_scan.tolist(), pad=0.015,
                  label=r"radiance [W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$]")
fig_scan.suptitle("Nadir radiance through the atmosphere at 800 km, thermal IR")

# The same scan as Planck brightness temperature per band.  The surface cells
# glow near their own 250-310 K, and the 1145 /cm panel, coupled to the ozone
# layer, sits systematically colder than the clean 895 /cm window.
tb_scan = np.stack([radiance_to_tb(scan[:, :, k], wn2hz(w))
                    for k, w in enumerate(BAND_WAVENUMBERS_CM)], axis=2)
fig_tb_scan, ax_tb_scan = plt.subplots(1, 2, figsize=(9.0, 4.0), sharex=True, sharey=True,
                                       constrained_layout=True)
for ax, k in zip(ax_tb_scan, range(len(BAND_WAVENUMBERS_CM))):
    mesh = ax.pcolormesh(LON_EDGES, LAT_EDGES, tb_scan[:, :, k], cmap="magma",
                         vmin=tb_scan.min(), vmax=tb_scan.max())
    ax.set_title(f"{BAND_WAVENUMBERS_CM[k]:.0f} /cm brightness temperature")
    ax.set_xlabel("Longitude [deg]")
    ax.set_xticks(np.arange(-180.0, 180.0, 45.0))
    ax.set_yticks(np.arange(-90.0, 90.0, 45.0))
    ax.set_box_aspect(0.5)
    mark_probes(ax)
ax_tb_scan[0].set_ylabel("Latitude [deg]")
fig_tb_scan.colorbar(mesh, ax=ax_tb_scan.tolist(), pad=0.015,
                     label="brightness temperature [K]")
fig_tb_scan.suptitle("Nadir brightness temperature through the atmosphere at 800 km, thermal IR")

fig_spectra, ax_spectra = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
for label, _, color, spectrum in PROBES:
    ax_spectra.plot(BAND_WAVENUMBERS_CM, spectrum, "o-", ms=7, lw=1.5, label=label, color=color)
    for w, v in zip(BAND_WAVENUMBERS_CM, spectrum):
        ax_spectra.annotate(f"{v:.2e}", (w, v), textcoords="offset points", xytext=(0, 10),
                            ha="center", fontsize=7, color=color)
ax_spectra.annotate("O$_3$ band between the bands", (1020.0, 1.2e-12), ha="center",
                    fontsize=7, color="0.45", rotation=90)
ax_spectra.set_xlabel("Wavenumber [1/cm]")
ax_spectra.set_ylabel("Radiance [W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$]")
ax_spectra.set_xlim(850.0, 1190.0)
ax_spectra.set_ylim(0.0, 5.0e-12)
ax_spectra.set_title("Red vs blue surface through the atmosphere, nadir at 800 km")
ax_spectra.grid(True, which="both", alpha=0.3)
ax_spectra.legend(framealpha=0.9, edgecolor="0.6")

fig_tb, ax_tb = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
for label, _, color, spectrum in PROBES:
    tb = radiance_to_tb(spectrum, wn2hz(BAND_WAVENUMBERS_CM))
    ax_tb.plot(BAND_WAVENUMBERS_CM, tb, "o-", ms=7, lw=1.5, label=label, color=color)
    for w, v in zip(BAND_WAVENUMBERS_CM, tb):
        ax_tb.annotate(f"{v:.1f} K", (w, v), textcoords="offset points", xytext=(0, 10),
                       ha="center", fontsize=7, color=color)
ax_tb.set_xlabel("Wavenumber [1/cm]")
ax_tb.set_ylabel("Brightness temperature [K]")
ax_tb.set_xlim(850.0, 1190.0)
ax_tb.set_ylim(250.0, 320.0)
ax_tb.set_title("Red vs blue surface through the atmosphere, nadir at 800 km")
ax_tb.grid(True, which="both", alpha=0.3)
ax_tb.legend(framealpha=0.9, edgecolor="0.6")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% References

ref_blue = [3.393952690350506e-12, 2.280326647232757e-12]
ref_red = [4.2966292192347805e-12, 2.489274689068086e-12]
ref_scan = [
    [
        [2.956539281431157e-12, 1.5873750304833577e-12],
        [2.278554643811394e-12, 1.4079441284107198e-12],
        [2.9565392814311566e-12, 1.5873750304833573e-12],
        [2.278554643811394e-12, 1.4079441284107194e-12],
        [2.9565392814311566e-12, 1.5873750304833575e-12],
        [2.2785546438113934e-12, 1.4079441284107194e-12],
        [2.956539281431157e-12, 1.5873750304833577e-12],
        [2.2785546438113942e-12, 1.4079441284107198e-12],
    ],
    [
        [3.393952690350506e-12, 2.280326647232757e-12],
        [4.2966292192347805e-12, 2.4892746890680863e-12],
        [3.393952690350506e-12, 2.2803266472327573e-12],
        [4.29662921923478e-12, 2.489274689068086e-12],
        [3.393952690350506e-12, 2.280326647232757e-12],
        [4.2966292192347805e-12, 2.489274689068086e-12],
        [3.3939526903505064e-12, 2.2803266472327577e-12],
        [4.2966292192347805e-12, 2.489274689068086e-12],
    ],
    [
        [4.2966292192347805e-12, 2.4892746890680887e-12],
        [3.393952690350506e-12, 2.280326647232757e-12],
        [4.2966292192347805e-12, 2.489274689068086e-12],
        [3.393952690350506e-12, 2.280326647232756e-12],
        [4.2966292192347805e-12, 2.489274689068085e-12],
        [3.393952690350506e-12, 2.2803266472327573e-12],
        [4.2966292192347805e-12, 2.489274689068086e-12],
        [3.3939526903505052e-12, 2.280326647232758e-12],
    ],
    [
        [2.2785546438113942e-12, 1.4079441284107198e-12],
        [2.9565392814311566e-12, 1.587375030483356e-12],
        [2.2785546438113934e-12, 1.4079441284107194e-12],
        [2.9565392814311566e-12, 1.5873750304833575e-12],
        [2.278554643811394e-12, 1.4079441284107194e-12],
        [2.9565392814311566e-12, 1.5873750304833567e-12],
        [2.2785546438113934e-12, 1.4079441284107194e-12],
        [2.956539281431157e-12, 1.5873750304833577e-12],
    ],
]

assert np.allclose(blue_spectrum, ref_blue, rtol=1e-8, atol=1e-30), \
    f"Mismatch to blue cell reference.\nReference:   {ref_blue}\nSimulation:  {blue_spectrum}"

assert np.allclose(red_spectrum, ref_red, rtol=1e-8, atol=1e-30), \
    f"Mismatch to red cell reference.\nReference:   {ref_red}\nSimulation:  {red_spectrum}"

assert np.allclose(scan, ref_scan, rtol=1e-8, atol=1e-30), \
    "Mismatch to the scan reference"

# The physics the plots show, pinned as inequalities so that a regression that
# breaks them fails even before the reference arrays are compared.  The red
# cells are warmer (+10 K) and emit better at 895 /cm, so they are brighter at
# both bands; at 1145 /cm the ozone band damps the surface contrast, so the
# red/blue ratio must be smaller there than at the clean 895 /cm window.  And
# the quadratic temperature profile must make the equatorial rows of the scan
# brighter than the polar ones -- each row averages 4 red and 4 blue cells, so
# the checkerboard and its temperature offset wash out of the row means.
for k in range(len(BAND_WAVENUMBERS_CM)):
    assert red_spectrum[k] > blue_spectrum[k], \
        f"Red not brighter than blue at {BAND_WAVENUMBERS_CM[k]} /cm"
assert red_spectrum[0] / blue_spectrum[0] > red_spectrum[1] / blue_spectrum[1], \
    "Red/blue contrast not damped at 1145 /cm by the ozone band"

row_means = scan.mean(axis=(1, 2))
assert row_means[1] > row_means[0], "22.5 N row not brighter than 67.5 N row"
assert row_means[2] > row_means[3], "22.5 S row not brighter than 67.5 S row"
