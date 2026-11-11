"""
Tests Profile2Path against the default path agenda, with and without wind.

With wind, Profile2Path shifts the frequencies of each path point by its
wind and interpolates the propagation matrices of the profile to them, so
radiances and Jacobians, those of the wind included, agree to within the
interpolation error.  Without wind, nothing is interpolated, so the
radiances are the same to the bit.  A profile computed without wind cannot
cover a path with wind, so Profile2Path throws unless f_extrapolfac allows
the extrapolation.

The propagation matrix and its temperature and frequency derivatives are
also compared directly with spectral_propmat_agenda, at one position and many
lines of sight.

Outside of ctest, that is without ARTS_HEADLESS, the differences are plotted.
"""

import os

import numpy as np
import pyarts3 as pyarts

PLOT = "ARTS_HEADLESS" not in os.environ
if PLOT:
    import matplotlib.pyplot as plt

    COLORS = ["#266CC1", "#eb6834", "#1baf7a"]  # Default and Profile2Path, or one per line of sight
    DEFAULT = {"color": COLORS[0], "label": "Default"}
    PROFILE = {"color": COLORS[1], "ls": "--", "label": "Profile2Path"}

    def angles(los):
        return f"{los[0]:g}°, {los[1]:g}°"

INF = float("inf")
POS = [600e3, 0.0, 0.0]
LOS = [[113.0, 30.0], [113.5, 60.0], [150.0, 200.0]]  # Limb at about 45 and 20 km, and down-looking
WIND = {"u": 100.0, "v": 50.0, "w": 5.0}

ws = pyarts.workspace.Workspace()

line_f0 = 118750348044.712
ws.freq_grid = np.linspace(-2e6, 2e6, 2001) + line_f0

ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=118e9, fmax=119e9)
ws.WignerInit()
ws.spectral_propmat_agendaAuto()

ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1)

ws.ray_path_observer_agendaSetGeometric(max_step_option="None", add_crossings=1, add_limb=1)


def set_wind(scale):
    for component, speed in WIND.items():
        ws.atm_field[f"wind_{component}"] = scale * speed


set_wind(1.0)
ws.measurement_sensors = []
ws.jac_targetsInit()
ws.jac_targetsAddTemperature()
ws.jac_targetsAddSpeciesVMR(species="O2")
for component in WIND:
    ws.jac_targetsAddWindField(component=component)
ws.jac_targetsFinalize()


def spectra():
    rad, jac = [], []
    for los in LOS:
        ws.spectral_rad_observer_agendaExecute(obs_pos=POS, obs_los=los)
        rad.append(np.array(ws.spectral_rad)[:, 0])
        jac.append(np.array(ws.spectral_rad_jac)[..., 0])
    return np.array(rad), np.array(jac)


def default_spectra():
    ws.spectral_propmat_and_atm_path_agendaSet(option="Default")
    return spectra()


def profile():
    """A profile at every altitude of the paths"""
    ws.spectral_propmat_and_atm_path_agendaSet(option="Default")
    alts = []
    for los in LOS:
        ws.ray_path_observer_agendaExecute(obs_pos=POS, obs_los=los)
        alts.extend(p.pos[0] for p in ws.ray_path)
    ws.alt_grid = np.unique(alts)
    ws.lat = 0.0
    ws.lon = 0.0
    ws.ProfileFromAltitude()
    ws.spectral_propmat_and_atm_path_agendaSet(option="Profile2Path")


# %% With wind

rad_ref, jac_ref = default_spectra()
profile()
rad, jac = spectra()

f = np.asarray(ws.freq_grid)
f_profile = np.asarray(ws.freq_grid_profile)
beta = np.linalg.norm(list(WIND.values())) / pyarts.arts.constants.c
assert f_profile[0] <= f[0] * (1 - beta) and f_profile[-1] >= f[-1] * (1 + beta), "The profile must cover the shift"
assert np.all(np.isin(f, f_profile)), "The profile must contain freq_grid"

# Errors relative to the largest value of each Jacobian element, so elements that no path reaches must agree at zero
rad_err = np.abs(rad - rad_ref).max() / np.abs(rad_ref).max()
jac_err = np.abs(jac - jac_ref).max(axis=(0, 2))
jac_max = np.abs(jac_ref).max(axis=(0, 2))
print(f"With wind: radiance error {rad_err:.1e}, largest relative Jacobian error {(jac_err[jac_max > 0] / jac_max[jac_max > 0]).max():.1e}")
print(f"  of the wind targets: {jac_err[-3:] / jac_max[-3:]}")

if PLOT:
    offset = (f - line_f0) / 1e6
    rj = pyarts.arts.constants.c**2 / (2 * pyarts.arts.constants.k * f**2)  # Radiance to Rayleigh-Jeans temperature

    fig, ax = plt.subplots(2, len(LOS), figsize=(13, 6), sharex=True, layout="constrained")
    for i, los in enumerate(LOS):
        ax[0, i].plot(offset, rad_ref[i] * rj, **DEFAULT)
        ax[0, i].plot(offset, rad[i] * rj, **PROFILE)
        ax[1, i].plot(offset, (rad[i] - rad_ref[i]) * rj, color=COLORS[1])
        ax[0, i].set_title(f"Line of sight {angles(los)}")
        ax[1, i].set_xlabel("Frequency - line centre [MHz]")
    ax[0, 0].set_ylabel("Brightness temperature (RJ) [K]")
    ax[1, 0].set_ylabel("Profile2Path - Default [K]")
    ax[0, 0].legend()
    fig.suptitle("Radiance with wind")
    plt.show()

    # The Jacobian of each target summed over its grid: the derivative by a change of the whole field
    targets = list(ws.jac_targets.atm)
    fig, ax = plt.subplots(len(targets), 2, figsize=(13, 2.6 * len(targets)), sharex=True, layout="constrained")
    for j, target in enumerate(targets):
        x = slice(target.x_start, target.x_start + target.x_size)
        for i, los in enumerate(LOS):
            ref_sum = jac_ref[i][x].sum(axis=0) * rj
            p2p_sum = jac[i][x].sum(axis=0) * rj
            ax[j, 0].plot(offset, ref_sum, color=COLORS[i], label=f"{angles(los)}, Default")
            ax[j, 0].plot(offset, p2p_sum, color=COLORS[i], ls="--", label=f"{angles(los)}, Profile2Path")
            ax[j, 1].plot(offset, p2p_sum - ref_sum, color=COLORS[i], label=angles(los))
        ax[j, 0].set_ylabel(f"{target.type}\n[K per unit]")
        ax[j, 1].set_ylabel("Profile2Path - Default")
    for a in ax[-1]:
        a.set_xlabel("Frequency - line centre [MHz]")
    ax[0, 0].legend(fontsize="small")
    ax[0, 1].legend(fontsize="small", title="Line of sight")
    fig.suptitle("Jacobians with wind, summed over the grid of each target")
    plt.show()

assert rad_err < 1e-4, rad_err
assert np.all(jac_err <= 1e-3 * jac_max), jac_err / jac_max

# %% The propagation matrix at one position, against spectral_propmat_agenda at many lines of sight
#
# The line of sight sets the wind shift, so it sets what Profile2Path interpolates.  The reference is
# spectral_propmat_agenda at the atmospheric point and shifted frequencies of the Default option.  The position is
# high, where the line is narrow, so the errors are large.  The wind rows of the Jacobian are the frequency derivative
# times the same factors for both, so their relative error is that of the frequency derivative.

alt_grid = np.asarray(ws.alt_grid)
POINT = [alt_grid[np.abs(alt_grid - 70e3).argmin()], 10.0, 20.0]  # Profile2Path needs an altitude of alt_grid
ANGLES = [[za, 30.0] for za in np.linspace(0, 180, 13)] + [[90.0, aa] for aa in np.linspace(0, 330, 12)]


def path_point(los):
    point = pyarts.arts.PropagationPathPoint()
    point.pos_type = pyarts.arts.PathPositionType("atm")
    point.los_type = pyarts.arts.PathPositionType("atm")
    point.pos = POINT
    point.los = los
    return point


def relative_error(x, ref):
    return np.abs(x - ref).max() / np.abs(ref).max()


ws.ray_path = [path_point(los) for los in ANGLES]
ws.Profile2Path()

atm_path = pyarts.arts.ArrayOfAtmPoint()
freq_grid_path = pyarts.arts.ArrayOfAscendingGrid()
freq_wind_shift_jac_path = pyarts.arts.ArrayOfVector3()
ws.atm_pathFromPath(atm_path=atm_path)
ws.freq_grid_pathFromPath(
    freq_grid_path=freq_grid_path, freq_wind_shift_jac_path=freq_wind_shift_jac_path, atm_path=atm_path
)

rows = {str(target.type): target.target_pos for target in ws.jac_targets.atm}
row_t = rows["t"]
rows_wind = [rows[f"wind_{component}"] for component in WIND]
df = f[1] - f[0]

print(f"Profile2Path against spectral_propmat_agenda at {POINT[0] / 1e3:g} km, relative to the largest reference value")
print(f"  {'zenith':>6} {'azimuth':>7} {'shift / df':>10} {'K':>8} {'dK/dT':>8} {'dK/df':>8}")
errors = []
shifts = []
curves = []  # Of the absorption (A) and its temperature and frequency derivatives, for the plots
for i, los in enumerate(ANGLES):
    assert np.array_equal(freq_grid_path[i], ws.freq_grid_path[i]), "Profile2Path must shift as the Default option"
    ws.spectral_propmat_agendaExecute(
        freq_grid=freq_grid_path[i],
        freq_wind_shift_jac=freq_wind_shift_jac_path[i],
        ray_point=ws.ray_path[i],
        atm_point=atm_path[i],
    )
    k_ref, dk_ref = np.array(ws.spectral_propmat), np.array(ws.spectral_propmat_jac)
    k, dk = np.array(ws.spectral_propmat_path[i]), np.array(ws.spectral_propmat_jac_path[i])

    errors.append(
        [relative_error(k, k_ref), relative_error(dk[row_t], dk_ref[row_t]), relative_error(dk[rows_wind], dk_ref[rows_wind])]
    )
    shifts.append((freq_grid_path[i][0] - f[0]) / df)
    print(f"  {los[0]:6.1f} {los[1]:7.1f} {shifts[-1]:10.3f} " + " ".join(f"{e:8.1e}" for e in errors[-1]))

    # The frequency derivative from the wind row of the component that shifts the most
    fs, shift_jac = np.asarray(freq_grid_path[i]), np.asarray(freq_wind_shift_jac_path[i])
    c = np.abs(shift_jac).argmax()
    to_df = 1 / (fs * shift_jac[c])
    curves.append(
        [
            [k_ref[:, 0], dk_ref[row_t, :, 0], dk_ref[rows_wind[c], :, 0] * to_df],
            [k[:, 0], dk[row_t, :, 0], dk[rows_wind[c], :, 0] * to_df],
        ]
    )

errors = np.array(errors)

if PLOT:
    names = ["K (A) [1/m]", "dK/dT [1/(m K)]", "dK/df [1/(m Hz)]"]

    # The shift is a whole number of grid spacings plus a fraction; the interpolation error grows with the fraction
    fraction = np.abs(np.array(shifts) - np.round(shifts))
    picks = [fraction.argmax(), np.abs(fraction - 0.25).argmin(), fraction.argmin()]

    fig, ax = plt.subplots(3, 2, figsize=(13, 9), sharex=True, layout="constrained")
    for q, name in enumerate(names):
        i = picks[0]
        offset = (np.asarray(freq_grid_path[i]) - line_f0) / 1e6
        ax[q, 0].plot(offset, curves[i][0][q], **DEFAULT)
        ax[q, 0].plot(offset, curves[i][1][q], **PROFILE)
        ax[q, 0].set_ylabel(name)
        for n, i in enumerate(picks):
            offset = (np.asarray(freq_grid_path[i]) - line_f0) / 1e6
            ax[q, 1].plot(
                offset,
                curves[i][1][q] - curves[i][0][q],
                color=COLORS[n],
                label=f"{angles(ANGLES[i])}, shift {shifts[i]:.3f} df",
            )
        ax[q, 1].set_ylabel("Profile2Path - Default")
    ax[0, 0].set_title(f"Line of sight {angles(ANGLES[picks[0]])}, shift {shifts[picks[0]]:.3f} df")
    ax[0, 0].legend()
    ax[0, 1].legend(fontsize="small", title="Line of sight")
    for a in ax[-1]:
        a.set_xlabel("Frequency - line centre [MHz]")
    fig.suptitle(f"Propagation matrix at {POINT[0] / 1e3:g} km")
    plt.show()

    fig, ax = plt.subplots(figsize=(7, 4.5), layout="constrained")
    for q, name in enumerate(names):
        ax.semilogy(fraction, errors[:, q], "o", ms=8, color=COLORS[q], label=name.split(" [")[0])
    ax.set_xlabel("Distance of the shift to a whole number of grid spacings [df]")
    ax.set_ylabel("Largest error / largest value")
    ax.legend()
    ax.set_title(f"Propagation matrix at {POINT[0] / 1e3:g} km, all {len(ANGLES)} lines of sight")
    plt.show()
assert np.all(errors < 1e-4), errors.max(axis=0)

# %% Without wind

set_wind(0.0)
rad_ref, jac_ref = default_spectra()
profile()
rad, jac = spectra()

assert np.array_equal(np.asarray(ws.freq_grid_profile), f), "Without wind, the profile is on freq_grid"
assert np.array_equal(rad, rad_ref), "Without wind, nothing is interpolated"
assert np.allclose(jac, jac_ref, rtol=1e-12, atol=0), "Without wind, nothing is interpolated"

# %% A path with wind through a profile without wind

set_wind(1.0)
try:
    spectra()
except RuntimeError:
    pass
else:
    raise AssertionError("Profile2Path must throw when the wind shifts frequencies beyond the profile")


@pyarts.workspace.arts_agenda(ws=ws, fix=True)
def spectral_propmat_and_atm_path_agenda(ws):
    ws.Profile2Path(f_extrapolfac=INF)
    ws.Touch(ws.ray_path)


rad, jac = spectra()
assert np.all(np.isfinite(rad)) and np.all(np.isfinite(jac)), "An infinite f_extrapolfac extrapolates"
