"""
What a checkerboard surface looks like from a high orbit through an atmosphere.

This is the atmospheric counterpart of ``1-surface-scattering-checkerboard.py``:
the same 4 x 8 checkerboard of red and blue Lambertian cells, the same
blackbody sun and the same Dirac sensor at 36 000 km nadir, but the empty
atmosphere is replaced by a tropical clearsky atmosphere up to 100 km.  The
air absorbs through H2O, O2 and CO2 lines and Rayleigh-scatters sunlight, so
the sensor now sees three contributions: sunlight reflected by the surface,
attenuated on its way up, and the sky itself glowing into the line of sight.

The absorption is not taken line-by-line but from an absorption lookup table
(*abs_lookup_data*): a precomputation of the line-by-line model on a grid of
altitude, temperature perturbation and humidity perturbation, which the ray
tracing then only interpolates.  The table is always recalculated at the start
of the script -- no caching -- because the example doubles as a regression
test and must be self-contained.  For the small grid used here (51 altitudes,
5 temperature and 5 humidity perturbations at 3 frequencies) that build takes
about a second, far less than the ray tracing itself.

With a sky to reflect, the surface agenda switches from the
SurfaceScatteringModelDirectOnly option of the atmosphere-free example to the
full SurfaceScatteringModel option: the diffuse sky arriving at a hit point is
integrated over a deliberately coarse quadrature hemisphere (3 zenith x 4
azimuth directions), while the direct solar beam is still counted exactly once
by the direct part of the chain.  The sky radiance along those quadrature
directions, and along the observer path itself, comes from first-order solar
Rayleigh scattering: the scattering propmat agenda AirSimple adds the air
scattering coefficient and the Rayleigh phase matrix to the absorption
propmat, and the incoming and observer agendas call
spectral_radClearskyRayleighScattering instead of the plain emission one.
Those two agendas are written out by hand rather than taken from the
predefined ClearskyRayleighScattering option, because that option computes
*ray_path_suns_path* -- the paths from every path point to the sun that the
first-order scattering source is integrated over -- only once for the observer
path.  The incoming agenda is re-executed for every quadrature direction at
every surface hit, so each of those rays needs its own sun paths recomputed
inside the agenda.

Two things are measured and plotted, as in the atmosphere-free example.  The
spectra above one blue and one red cell show the reflectivity spectra
imprinted on the upwelling radiance, and how the band with r = 0 no longer
collapses to the thermal floor but to the brightness of the sky.  The scan
over the centre of every cell shows the checkerboard modulated by the day-
night pattern of the sun; first-order scattering needs sunlight along the
path, so the night side stays dark at every band.

The example ends by pinning the spectra and the scan to reference arrays and
by checking that exactly the sunlit cells are bright, so it doubles as a
regression test of the framework.
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

# Ascending in frequency means descending in wavelength
BAND_WAVELENGTHS_NM = [650.0, 550.0, 400.0]
REFLECTIVITY_WAVELENGTHS_NM = [700.0, 600.0, 500.0, 400.0]
RED_REFLECTIVITY = [1.0, 0.9, 0.5, 0.0]
BLUE_REFLECTIVITY = [0.0, 0.25, 0.7, 1.0]

CELL_SIZE = 45.0
SURF_T = 280.0
TOA = 100e3
MAX_STEPSIZE = 9600.0

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

# Radiance floors of the plots: values at or below are not resolved.  The
# spectra floor sits just below the airlight of the r = 0 band, so that the
# sky glow shows up as the new dark limit instead of the thermal one.
RADIANCE_FLOOR = 1e-11
SPECTRA_FLOOR = 1e-15

# Absorbing species and the line-selection window covering the three
# measurement bands (461.2, 545.0, and 749.5 THz) with a ~1 THz margin on
# each side
ABS_SPECIES = ["H2O-161", "O2-66", "CO2-626"]
BAND_FMIN = 460e12
BAND_FMAX = 751e12

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

    ws.freq_grid = nm2hz(BAND_WAVELENGTHS_NM)

    # Absorbing species and the line bands covering the three measurement bands
    ws.abs_speciesSet(species=ABS_SPECIES)
    ws.ReadCatalogData()
    ws.abs_bandsSelectFrequencyByLine(fmin=BAND_FMIN, fmax=BAND_FMAX)

    # Tropical clearsky atmosphere up to TOA (replaces the empty atm_fieldInit)
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


# Floor the radiance at RADIANCE_FLOOR so that the invisible night side does
# not stretch the log colour scale decades beyond the day-side glow
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
fig_scan.suptitle("Nadir radiance through the atmosphere at 36 000 km")

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
ax_spectra.annotate("r = 0 → airlight", (400.0, SPECTRA_FLOOR * 1.5), ha="center",
                    fontsize=7, color="0.45")
ax_spectra.set_xlabel("Wavelength [nm]")
ax_spectra.set_ylabel("Radiance [W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$]")
ax_spectra.set_ylim(SPECTRA_FLOOR, 3e-8)
ax_spectra.set_title("Red vs blue surface through the atmosphere, nadir at 36 000 km")
ax_spectra.grid(True, which="both", alpha=0.3)
ax_spectra.legend(framealpha=0.9, edgecolor="0.6")

if "ARTS_HEADLESS" not in os.environ:
    plt.show()

# %% References

ref_blue = [2.687358754136667e-09, 7.768156084590544e-09, 4.320857776343119e-09]
ref_red = [5.883894616607171e-09, 2.874498668959128e-09, 1.6858576347319554e-14]
ref_scan = [
    [
        [6.769403577754206e-18, 9.199325467754575e-18, 8.699342615530867e-18],
        [2.372526655996423e-10, 4.906866538188595e-10, 3.927524793112714e-11],
        [9.760994727462792e-09, 5.070912302369151e-09, 1.825812074754026e-14],
        [1.727776142628184e-09, 4.8576727137270195e-09, 2.3135367080559025e-09],
        [1.2695147825117281e-08, 6.758525209899363e-09, 1.9706843937299886e-14],
        [7.186548767212475e-10, 1.8394579457247626e-09, 5.10415801729797e-10],
        [1.1189078632215253e-14, 1.4918947444597076e-14, 1.2968584829906472e-14],
        [3.468195116959678e-23, 2.199169453476818e-22, 4.1132092338933274e-22],
    ],
    [
        [2.09312395917317e-40, 9.088114173113566e-46, 5.693125165918009e-59],
        [1.0314807922285136e-16, 1.3517454634022003e-16, 1.0831246439198079e-16],
        [1.8297365436863031e-09, 5.166191048973635e-09, 2.5198430671647926e-09],
        [2.6058701451730614e-08, 1.4510909944102396e-08, 3.058818673188402e-14],
        [2.687358754136667e-09, 7.768156084590544e-09, 4.320857776343119e-09],
        [5.883894616607171e-09, 2.874498668959128e-09, 1.6858576347319554e-14],
        [2.0931239591731086e-40, 9.088114173113213e-46, 5.693125165917806e-59],
        [3.3487216663924825e-40, 1.1102281625686992e-45, 2.3129612350352074e-59],
    ],
    [
        [3.3487216663929262e-40, 1.1102281625688821e-45, 2.3129612350360927e-59],
        [2.1291000728214465e-35, 3.031745798841404e-35, 9.12557475286055e-39],
        [1.0893754825402953e-08, 5.7199959064962965e-09, 1.8777595447416167e-14],
        [2.6393710735538233e-09, 7.622586910014014e-09, 4.218752763339413e-09],
        [1.7986200708806313e-08, 9.821742855246891e-09, 2.317197816258711e-14],
        [2.2229963784974836e-10, 4.5259585011385013e-10, 3.304815691378165e-11],
        [3.3487216663925763e-40, 1.1102281625687458e-45, 2.3129612350354944e-59],
        [2.0931239591736115e-40, 9.088114173115393e-46, 5.69312516591891e-59],
    ],
    [
        [2.093123959173447e-40, 9.088114173114725e-46, 5.693125165918523e-59],
        [1.7003784846365063e-34, 4.7729312761389304e-35, 2.3129612350357323e-59],
        [3.268727557450812e-11, 3.9837717018849354e-11, 1.6131898021011187e-13],
        [4.11873963136287e-09, 1.9007644508644713e-09, 1.6411442959598655e-14],
        [3.1515083643578675e-10, 6.95252585525419e-10, 8.125786935859717e-11],
        [8.135573589808023e-17, 1.1345357822847836e-16, 1.0786019946714246e-16],
        [2.0931239591734532e-40, 9.088114173114817e-46, 5.693125165918596e-59],
        [3.3487216663927656e-40, 1.1102281625688162e-45, 2.3129612350357345e-59],
    ],
]

# The sub-floor thermal terms are only loosely pinned (atol dwarfs them); the
# reflected sunlight and the airlight, which is what the plots show, are
# pinned to relative machine precision
assert np.allclose(blue_spectrum, ref_blue, rtol=1e-8, atol=1e-30), \
    f"Mismatch to blue cell reference.\nReference:   {ref_blue}\nSimulation:  {blue_spectrum}"

assert np.allclose(red_spectrum, ref_red, rtol=1e-8, atol=1e-30), \
    f"Mismatch to red cell reference.\nReference:   {ref_red}\nSimulation:  {red_spectrum}"

assert np.allclose(scan, ref_scan, rtol=1e-8, atol=1e-30), \
    "Mismatch to the scan reference"

# The scan is the checkerboard modulated by the sun: every sunlit cell glows
# with reflected sunlight and airlight above the floor, every night-side cell
# stays below it -- first-order Rayleigh scattering needs sunlight along the
# path, and the thermal emission of the atmosphere is nothing at these
# frequencies.  The two cells closest to the terminator sit just on either
# side of it ((67.5, 112.5) on the night side at 1e-14, (-67.5, -67.5) on the
# day side at 1e-11), and the 1e-12 tolerance on the cosine of the incidence
# angle keeps them unambiguously on the side their radiance puts them on
for i, lat in enumerate(CELL_LATS):
    for j, lon in enumerate(CELL_LONS):
        cos_inc = (np.sin(np.deg2rad(SUN_LAT)) * np.sin(np.deg2rad(lat)) +
                   np.cos(np.deg2rad(SUN_LAT)) * np.cos(np.deg2rad(lat)) *
                   np.cos(np.deg2rad(lon - SUN_LON)))
        sunlit = cos_inc > 1e-12
        bright = bool(np.any(scan[i, j] > RADIANCE_FLOOR))
        assert bright == sunlit, f"Cell ({lat}, {lon}) sunlit={sunlit}, scan={scan[i, j]}"
