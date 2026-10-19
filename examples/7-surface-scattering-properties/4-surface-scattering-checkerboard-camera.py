"""
A camera photograph of the checkerboard surface through the atmosphere.

This is the camera counterpart of ``2-surface-scattering-checkerboard-atm-lut.py``:
the same 4 x 8 checkerboard of red and blue Lambertian cells, the same tropical
clearsky atmosphere up to 100 km with the absorption taken from a lookup table,
and the same first-order solar Rayleigh scattering lighting the sky.  What
changes is the instrument: instead of pointing a sensor at one spot after
another, a thin-lens camera with a 25 x 13 pixel CCD looks nadir at the surface
from 36 000 km, and one image comes out of a single simulation.

The measurement is taken at a single band, 550 nm, unlike the three bands of
the previous examples.  One frequency keeps the run short enough for the test
suite -- the camera traces one full observer chain per pixel -- and it is also
all that is needed: with one band there is no colour to render, the image is
the radiance itself, plotted as a single radiance map.  The band sits where
the two surfaces still contrast (red r = 0.9, blue r = 0.25), and the
reflectivity spectra carry a node exactly at 550 nm so the band samples those
values without interpolation.

The camera is described by its CCD (16 mm x 9 mm), a lens (27 mm focal length,
0.1 m aperture) and a focus distance (the surface, at infinity for the lens
equation).  Each pixel maps to one line of sight through the thin-lens formula,
so the field of view follows from the CCD and the focal length: the half-angle
of about 16 deg in the vertical direction reaches past the geometric horizon
of the Earth seen from 36 000 km, which sits at asin(R/(R + h)) ~ 8.6 deg.
Inside that cone the rays graze the surface and the ground point runs off
towards the limb, so the image shows the whole disk of the checkerboard,
several cells across, framed by the black of space.  The camera looks
nadir at (0 N, 45 E) -- exactly a corner of the checkerboard -- and the sub-
solar point sits at (11.25 N, 11.25 W) as in the previous examples, so the
image carries the full day-night pattern: the cells near the subsolar point
glow, the far side of the disk fades into the night.

The image comes back in the flat *measurement_vec*, one value per pixel and
frequency in row-major order; *measurement_sensor_metaFromMeasurementVec*
reshapes it into a gridded field whose grids are the angular offsets of the
pixels from the central line of sight.  Those grids are what the image plot
uses for its axes.  As a cross-check, the scene is also measured with a plain
nadir sensor at the camera's position: with an odd pixel count the central
pixel of the camera sits exactly on the camera's line of sight, so that
sensor follows the same ray and must reproduce the central pixel exactly.

The example ends by pinning the image and the centre measurement to reference
arrays and by checking the scene physically -- the corner pixels look past the
horizon into space, the centre pixel sees the sunlit cell corner, the disk
carries the red/blue contrast, and the side of the disk turned towards the
sun is brighter than the far side -- so it doubles as a regression test of the
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

# A single band at 550 nm; the reflectivity spectra carry a node here, so the
# red and blue cells sample exactly 0.9 and 0.25
BAND_WAVELENGTH_NM = 550.0
REFLECTIVITY_WAVELENGTHS_NM = [700.0, 600.0, 550.0, 500.0, 400.0]
RED_REFLECTIVITY = [1.0, 0.9, 0.9, 0.5, 0.0]
BLUE_REFLECTIVITY = [0.0, 0.25, 0.25, 0.7, 1.0]

CELL_SIZE = 45.0
SURF_T = 280.0
TOA = 100e3
MAX_STEPSIZE = 9600.0

SUN_LAT, SUN_LON = 11.25, -11.25
SENSOR_ALT = 36e6
SENSOR_LAT, SENSOR_LON = 0.0, 45.0

# The geometric horizon of the Earth seen from the sensor altitude: rays inside
# this angular offset from nadir hit the surface, rays outside leave it for
# space.  The camera's 16.5 deg half-FOV reaches past it, so the frame corners
# are black.
EARTH_RADIUS = 6371e3
HORIZON_OFFSET = np.rad2deg(np.arcsin(EARTH_RADIUS / (EARTH_RADIUS + SENSOR_ALT)))

# The camera: 25 x 13 pixels on a 16 mm x 9 mm CCD behind a 27 mm lens.  The
# odd pixel counts put one pixel exactly on the central line of sight, which
# makes the nadir cross-check below exact.  The vertical half-FOV is
# atan(8/27) = 16.5 deg, past the ~8.6 deg geometric horizon of the Earth from
# 36 000 km, so the frame holds the full disk plus the black of space around
# it.  The pixel count is deliberately small: every pixel on the disk costs a
# full observer chain with the quadrature hemisphere at its surface hit, and
# the example doubles as a test run.  The focus distance is the surface
# distance, so the whole scene is in focus.
N_H, N_W = 25, 13
CCD_H, CCD_W = 0.016, 0.009
FOCAL_LENGTH = 0.027
APERTURE_DIAMETER = 0.1
FOCUS_DISTANCE = SENSOR_ALT

# Radiance floor of the log image plot: values at or below are not resolved.
# It sits below the sunlit surface glow but well above the thermal floor, so
# the airlight and the dark cells still show structure.
RADIANCE_FLOOR = 1e-12

# Absorbing species and the line-selection window covering the measurement
# band (545.1 THz) with a ~1 THz margin on each side
ABS_SPECIES = ["H2O-161", "O2-66", "CO2-626"]
BAND_FMIN = 544e12
BAND_FMAX = 546e12

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


def nm2hz(wavelength_nm):
    """Wavelength [nm] to frequency [Hz]"""
    return SPEED_OF_LIGHT / (np.asarray(wavelength_nm, dtype=float) * 1e-9)


# %% The two surface scattering models

def reflectivity_spectrum(values):
    """Reflectivity spectrum on its own ascending frequency grid"""
    return arts.SpectralGeodeticField3(
        name="reflectivity",
        grid_names=["Frequency", "Latitude", "Longitude"],
        grids=(nm2hz(REFLECTIVITY_WAVELENGTHS_NM), [0.0], [0.0]),
        data=np.asarray(values, dtype=float).reshape(-1, 1, 1),
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

    ws.freq_grid = np.array([nm2hz(BAND_WAVELENGTH_NM)])

    # Absorbing species and the line band covering the measurement band
    ws.abs_speciesSet(species=ABS_SPECIES)
    ws.ReadCatalogData()
    ws.abs_bandsSelectFrequencyByLine(fmin=BAND_FMIN, fmax=BAND_FMAX)

    # Tropical clearsky atmosphere up to TOA (replaces the empty atm_fieldInit)
    ws.atm_fieldRead(toa=TOA, basename="planets/Earth/afgl/tropical/", missing_is_zero=1)

    # The absorption lookup table, always recalculated here at the start of the
    # script so that the example carries no cached state.  This precomputes the
    # line-by-line absorption on the altitude x temperature x humidity grid of
    # the constants above; the ray tracing below only interpolates it.  With a
    # single frequency the build is a fraction of a second.
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

    # The Rayleigh methods require an (empty) jac_targets
    ws.jac_targets = arts.JacobianTargets()

    # Rayleigh scattering propmat for air, added on top of the absorption one
    ws.spectral_propmat_scat_agendaSet(option="AirSimple")

    ws.surf_fieldEarth()
    ws.surf_field["t"] = SURF_T
    ws.surf_field[RED] = red_mask
    ws.surf_field[BLUE] = lambda lat, lon: 1.0 - red_mask(lat, lon)

    ws.surface_models = surface_models()

    ws.sunBlackbody(lat=SUN_LAT, lon=SUN_LON)
    ws.suns = [ws.sun]

    ws.zen_grid, ws.az_grid, ws.zen_grid_weights, ws.az_grid_weights = quadrature_grids()

    ws.max_stepsize = MAX_STEPSIZE

    ws.ray_path_observer_agendaSetGeometric()
    ws.spectral_rad_space_agendaSet(option="SunOrCosmicBackground")
    ws.spectral_rad_closed_surface_agendaSet(option="Blackbody")
    # With an atmosphere there is a sky to reflect, so the surface agenda takes
    # the full option: the diffuse sky is integrated over the quadrature
    # hemisphere, and the direct solar beam is added exactly once by the direct
    # part of the chain (quadrature directions containing a sun are dropped).
    ws.spectral_rad_surface_agendaSet(option="SurfaceScatteringModel")

    # Custom incoming/observer agendas with solar Rayleigh scattering.
    #
    # The predefined spectral_rad_incoming_agendaSet(option="ClearskyRayleighScattering")
    # never recomputes ray_path_suns_pathFromPathObserver per ray, so the sun
    # path would be stale for every ray after the first;
    # spectral_radClearskyRayleighScattering consumes ray_path_suns_path as a
    # hard input.  Recomputing the sun path inside the agendas keeps it
    # consistent with the current path.
    @pa.workspace.arts_agenda(ws=ws, fix=True)
    def spectral_rad_incoming_agenda(ws):     # downwelling radiance at surface hits
        ws.ray_path_observer_agendaExecute()
        ws.ray_path_suns_pathFromPathObserver(just_hit=1)
        ws.spectral_radClearskyRayleighScattering()

    @pa.workspace.arts_agenda(ws=ws, fix=True)
    def spectral_rad_observer_agenda(ws):     # sensor path incl. airlight
        ws.ray_path_observer_agendaExecute()
        ws.ray_path_suns_pathFromPathObserver(just_hit=1)
        ws.spectral_radClearskyRayleighScattering()
        ws.spectral_rad_jacAddSensorJacobianPerturbations()

    return ws


# %% The camera image: one exposure of the checkerboard at 550 nm

ws = build_workspace()

ws.measurement_sensorInit()
ws.measurement_sensorAddCamera(pos=[SENSOR_ALT, SENSOR_LAT, SENSOR_LON],
                               los=[180.0, 0.0],
                               n_h=N_H,
                               n_w=N_W,
                               ccd_h=CCD_H,
                               ccd_w=CCD_W,
                               focal_length=FOCAL_LENGTH,
                               aperture_diameter=APERTURE_DIAMETER,
                               focus_distance=FOCUS_DISTANCE)
ws.measurement_vecFromSensor(kernel="High Performance")

# The flat measurement vector is one value per pixel and frequency, row-major
# over the CCD; the meta gridded fields hold the same values with the angular
# offsets of the pixels from the central line of sight as their grids.
ws.measurement_sensor_metaFromMeasurementVec()
image_gf = ws.measurement_sensor_meta[0].data
image = np.array(image_gf.data)[:, :, 0]
row_offsets = np.array(image_gf.grids[0])
col_offsets = np.array(image_gf.grids[1])

assert image.shape == (N_H, N_W)
assert np.allclose(np.array(ws.measurement_vec).reshape(N_H, N_W), image)

# %% The centre check: a plain nadir sensor on the central pixel must agree

ws.measurement_sensorSimple(pos=[SENSOR_ALT, SENSOR_LAT, SENSOR_LON], los=[180.0, 0.0])
ws.measurement_vecFromSensor(kernel="High Performance")
centre_spectrum = np.array(ws.measurement_vec)

# The central line of sight of the camera hits (0 N, 45 E) straight down, the
# same ray the nadir sensor traces, and with the odd pixel counts the central
# pixel sits exactly on that line of sight, so the two must agree to relative
# machine precision
centre_pixel = image[N_H // 2, N_W // 2]
assert np.allclose(centre_spectrum, centre_pixel, rtol=1e-8, atol=1e-30), \
    f"Camera centre pixel and nadir sensor disagree.\nCentre pixel:   {centre_pixel}\nNadir sensor:   {centre_spectrum}"


# %% Plots

# The image on a log scale floored at RADIANCE_FLOOR: the sunlit cells glow
# within a factor of a few of each other, while the night side and the space
# outside the horizon sit many decades below
log_image = np.log10(np.maximum(image, RADIANCE_FLOOR))
fig_image, ax_image = plt.subplots(figsize=(7.0, 6.0), constrained_layout=True)
mesh = ax_image.pcolormesh(col_offsets, row_offsets, log_image, cmap="magma",
                           vmin=np.log10(RADIANCE_FLOOR), vmax=log_image.max())
for edge in HORIZON_OFFSET, -HORIZON_OFFSET:
    ax_image.axhline(edge, color="cyan", ls="--", lw=0.8)
    ax_image.axvline(edge, color="cyan", ls="--", lw=0.8)
ax_image.annotate("geometric horizon", (HORIZON_OFFSET, -HORIZON_OFFSET),
                  textcoords="offset points", xytext=(4, 6), fontsize=8, color="cyan")
ax_image.set_xlabel("Column angular offset [deg]")
ax_image.set_ylabel("Row angular offset [deg]")
ax_image.set_title("Checkerboard at 550 nm, camera at 36 000 km nadir")
ax_image.set_aspect("equal")
fig_image.colorbar(mesh, ax=ax_image, pad=0.02,
                   label=r"log$_{10}$ radiance [W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$]")

# A cut through the middle row of the image, against the column angular
# offset: the checkerboard cells show up as steps in radiance, the day-night
# pattern as a slow envelope, and the drop past the horizon offset as the
# corner rays leaving the Earth for space
cut_row = N_H // 2
fig_cut, ax_cut = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
ax_cut.semilogy(col_offsets, np.maximum(image[cut_row], 1e-20), "o-", ms=3, lw=1.2,
                color="tab:blue")
ax_cut.axvline(HORIZON_OFFSET, color="0.5", ls="--", lw=0.8)
ax_cut.axvline(-HORIZON_OFFSET, color="0.5", ls="--", lw=0.8)
ax_cut.set_xlabel("Column angular offset [deg]")
ax_cut.set_ylabel("Radiance [W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$]")
ax_cut.set_title(f"Row {cut_row} of the camera image, dashed lines at the horizon")
ax_cut.grid(True, which="both", alpha=0.3)

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% References

ref_centre = [9.54733602996089e-09]
ref_image = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 1.3814522714756874e-10, 9.577728564933682e-10, 1.6954022707370102e-09, 7.721312335333699e-09, 1.4668338109408261e-09, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 1.6161768479912645e-11, 8.158413055092572e-10, 5.9940122561544764e-09, 8.626622209056705e-09, 2.996296982402951e-09, 1.2338024235804459e-08, 1.2544536988568101e-08, 0.0, 0.0, 0.0],
    [0.0, 0.0, 1.4891966123707984e-16, 9.593248515216956e-10, 4.099295856627465e-09, 7.042997167174593e-09, 9.638609766720382e-09, 3.3074391253432403e-09, 3.819594133328575e-09, 4.1639659835804855e-09, 1.4716349565657828e-08, 0.0, 0.0],
    [0.0, 1.0269596869258372e-45, 1.2218366182846161e-14, 1.534982239824819e-09, 4.625794047730711e-09, 7.487165802802939e-09, 1.0074272632408689e-08, 3.429625668871837e-09, 3.978096473920051e-09, 4.399319153360705e-09, 4.6032059151673006e-09, 1.3952929155409476e-08, 0.0],
    [0.0, 1.882737280870815e-36, 1.822228856460961e-14, 1.7291066839197115e-09, 4.762217197006012e-09, 7.591983086558511e-09, 1.0145931501806423e-08, 3.4590150115735753e-09, 4.015488379727509e-09, 4.466549610593967e-09, 4.767633283151262e-09, 1.666420682982089e-08, 0.0],
    [0.0, 6.294005499536599e-36, 1.798694661561451e-14, 1.632272970979925e-09, 4.607902860024619e-09, 7.412014391655524e-09, 9.969657204804683e-09, 3.407760332066955e-09, 3.9694097679269715e-09, 4.439587768191659e-09, 4.771025557895947e-09, 1.7209696569801564e-08, 0.0],
    [1.077779397661542e-58, 2.011325604427297e-36, 8.56024158359916e-15, 1.2937072682277028e-09, 4.20176053201189e-09, 6.996014444486126e-09, 9.54733602996089e-09, 3.2911722089700754e-09, 3.850182240924115e-09, 4.321308363833004e-09, 4.663987282678399e-09, 1.6971336886296635e-08, 1.0016886844794007e-16],
    [0.0, 1.7323653510801796e-37, 1.5278928465892512e-15, 2.1556859713931347e-10, 9.919537918002203e-10, 1.7635761798376194e-09, 2.4716946680385857e-09, 1.1189405538264014e-08, 1.3212986424323734e-08, 1.4901875002206537e-08, 1.609138590786084e-08, 4.487627900765331e-09, 0.0],
    [0.0, 1.2548735619536021e-45, 8.180810144671929e-18, 5.8044912878632235e-11, 7.517223669464333e-10, 1.5183134567760714e-09, 2.2236227982288966e-09, 1.0294194019404006e-08, 1.2302221368989715e-08, 1.3924089370568561e-08, 1.5005570949981594e-08, 4.054808823530926e-09, 0.0],
    [0.0, 1.4098595595125885e-45, 8.677541723552777e-35, 1.9098376565810873e-14, 4.486250282049768e-10, 1.198357198050893e-09, 1.9050584317541054e-09, 9.115629427170413e-09, 1.1086271888925792e-08, 1.2610875483940159e-08, 1.3381600341670087e-08, 3.108591347219882e-09, 0.0],
    [0.0, 0.0, 3.7812847877336907e-39, 2.164526773815445e-16, 1.1401300835282376e-10, 7.936165269119163e-10, 1.490831272873526e-09, 7.602276994280658e-09, 9.449296072630042e-09, 1.073457748798541e-08, 2.965821484835278e-09, 0.0, 0.0],
    [0.0, 0.0, 0.0, 3.6771602632050524e-35, 1.8063766411397192e-14, 2.865496731685732e-10, 9.39127363394375e-10, 5.492414847924833e-09, 1.9684518387067235e-09, 2.101702021309969e-09, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 2.137344636322058e-37, 1.2272902499216133e-14, 6.187889985553975e-10, 5.859274602235607e-10, 1.645931402143617e-09, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
]

assert np.allclose(centre_spectrum, ref_centre, rtol=1e-8, atol=1e-30), \
    f"Mismatch to centre reference.\nReference:   {ref_centre}\nSimulation:  {centre_spectrum}"

assert np.allclose(image, ref_image, rtol=1e-8, atol=1e-30), \
    "Mismatch to the image reference"

# The physics the plots show, pinned as inequalities so that a regression that
# breaks them fails even before the reference arrays are compared.  The four
# corners of the frame look past the horizon into the cosmic background and
# stay dark.  The centre pixel sees the sunlit cell corner at (0 N, 45 E) and
# glows.  Inside the horizon cone the disk holds the full checkerboard, so its
# brightest and darkest surface pixels differ by the reflectivity and sun-
# angle contrast.  And the subsolar point at (11.25 N, 11.25 W) sits west of
# the nadir (0 N, 45 E), so the half of the disk turned towards the sun is
# clearly brighter than the far half.
corners = np.array([image[0, 0], image[0, -1], image[-1, 0], image[-1, -1]])
assert np.all(corners < 1e-15), f"Corner pixels not space-dark: {corners}"

assert centre_pixel > 1e-10, f"Centre pixel too dark: {centre_pixel}"

# The disk pixels are those whose total angular offset from nadir stays inside
# the horizon cone -- a circular mask, the square corners of the frame are
# space.  Inside it the checkerboard and the day-night pattern span more than
# two decades of radiance.
offset_r = np.sqrt(row_offsets[:, None] ** 2 + col_offsets[None, :] ** 2)
on_disk = offset_r < HORIZON_OFFSET
disk = image[on_disk]
assert disk.max() / disk.min() > 3.0, \
    f"Disk radiance contrast too small: {disk.max():.3e} / {disk.min():.3e}"

half_a = image[:, :N_W // 2][on_disk[:, :N_W // 2]].mean()
half_b = image[:, N_W // 2:][on_disk[:, N_W // 2:]].mean()
assert max(half_a, half_b) / min(half_a, half_b) > 1.5, \
    f"Day-night asymmetry too small: {half_a:.3e} vs {half_b:.3e}"
