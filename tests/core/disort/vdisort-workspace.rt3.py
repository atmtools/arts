"""The VDISORT workspace methods against references derived outside VDISORT.

vdisort_settings_agendaSetup and vdisort_settings_agendaExecute build
VDISORT's input from ARTS data, and vdisort_spectral_rad_fieldCalc maps its
output onto the radiance field.  Three references check them:

1. Clear sky: the exact solution of the transfer equation through layers whose
   Planck function is linear in optical depth, recomputed here from the path's
   own propagation matrices and temperatures.  The sky is given a polarization
   to check Q, U and V; the blackbody surface keeps them zero upward.  The top
   layer has no absorption, so min_optical_depth sets its optical thickness.
   The solver writes the source as c0 + c1 tau in the global optical depth tau,
   so the tolerance of I is 16 eps times its largest term, |c0| + |c1| (tau + 1),
   or the largest boundary value.  Q, U and V are the sky's, transmitted layer
   by layer, so theirs is 16 eps per layer times the sky's value.

2. Single scattering of the sun by an optically thin Henyey-Greenstein layer:
   the first-order closed form in the scattering angle between the beam,
   which propagates away from the sun, and the propagation direction (zenith,
   azimuth) of the field.  This pins the field's grids to propagation
   directions and the beam to the sun's position.  Higher orders are below
   q / (1 - q) of the largest first-order radiance, with
   q = omega tau mean|P| / mu_min.

3. Evans' doubling-adding RT3 (polradtran) on the same path, through
   rt3.problem_from_path: thermal emission with Rayleigh gas scattering and
   Mie drops above a blackbody surface, and the sun with Rayleigh, Henyey-
   Greenstein and Mie scattering above a Lambertian surface, in 8 Fourier
   modes and all four Stokes components.  RT3 differs by its first-order
   doubling error, below max_delta_tau / mu0 of max I, and must agree to 10
   times that.

The sign of U and V relative to ARTS's sensor polarization is not tested.
"""

import numpy as np
import pyarts3 as pa

A = pa.arts

FREQ = [89e9, 150e9]
NMU = 8  # streams per hemisphere
LEGENDRE = A.rt3.max_legendre_degree(NMU, A.polradtran.QuadratureType.double_gauss) + 1
MAX_DELTA_TAU = 1e-7
T_SURF = 288.0
EPS = np.finfo(float).eps

# Levels from the top; the gases vanish above 6500 m, which leaves the top layer empty
Z = np.array([8000.0, 7000.0, 6000.0, 4500.0, 3200.0, 2000.0, 1000.0, 0.0])
GAS = Z < 6500.0

hg_ext = A.ScatteringSpeciesProperty("hg", A.ParticulateProperty.Extinction)
hg_ssa = A.ScatteringSpeciesProperty("hg", A.ParticulateProperty.SingleScatteringAlbedo)
drops = A.ScatteringSpeciesProperty("drops", A.ParticulateProperty.NumberDensity)

ATMOSPHERE = {
    "t": 290.0 - 0.01 * Z,
    "p": 1e5 * np.exp(-Z / 8000.0),
    "H2O": np.where(GAS, 0.02 * np.exp(-Z / 2000.0), 0.0),
    "O2": np.where(GAS, 0.21, 0.0),
    hg_ext: 2e-4 * np.exp(-(((Z - 3000.0) / 1000.0) ** 2)),
    hg_ssa: np.full_like(Z, 0.9),
    drops: 30.0 * np.exp(-(((Z - 2500.0) / 1500.0) ** 2)),
}


def profile(values):
    field = A.GeodeticField3(
        grids=[Z[::-1], [0.0], [0.0]],
        grid_names=["altitude", "latitude", "longitude"],
        data=np.asarray(values, dtype=float)[::-1].reshape(len(Z), 1, 1),
    )
    return field


def workspace(species, gas=True, fields={}):
    """A workspace with the atmosphere (fields replace its profiles) and a path down through its levels."""
    ws = pa.Workspace()
    ws.freq_grid = FREQ
    ws.abs_speciesSet(species=["H2O-PWR98", "O2-PWR98"] if gas else [])
    if gas:
        ws.ReadCatalogData()
    ws.spectral_propmat_agendaAuto()
    ws.surf_fieldPlanet(option="Earth")
    ws.surf_field[A.SurfaceKey("t")] = T_SURF
    ws.atm_field.top_of_atmosphere = Z[0]
    for key, values in {**ATMOSPHERE, **fields}.items():
        if isinstance(key, str) and key in ("H2O", "O2"):
            key = A.SpeciesEnum(key)
        ws.atm_field[key] = profile(values)
        ws.atm_field[key].alt_low = "Linear"
        ws.atm_field[key].alt_upp = "Linear"
    ws.scat_species = species
    ws.spectral_propmat_scat_spectral_agendaSet(option="FromSpeciesTRO")

    path = A.ArrayOfPropagationPathPoint()
    for z in Z:
        point = A.PropagationPathPoint()
        point.pos = [z, 0.0, 0.0]
        point.los = [0.0, 180.0]
        point.pos_type = point.los_type = A.PathPositionType("atm")
        path.append(point)
    ws.ray_path = path
    ws.quadrature_dimension = 2 * NMU
    return ws


def field(ws, azimuths=(0.0,)):
    """Solve; return the field's data [freq, layer, azimuth, zenith, stokes] and propagation cosines."""
    ws.vdisort_settings_agendaExecute()
    ws.vdisort_spectral_rad_fieldCalc(azi_grid=list(azimuths))
    f = ws.vdisort_spectral_rad_field
    return np.asarray(f.data), np.cos(np.radians(np.asarray(f.zen_grid))), np.asarray(f.azi_grid)


def planck(f, t):
    h, k, c = A.constants.h, A.constants.k, A.constants.c
    return 2 * h * f**3 / c**2 / np.expm1(h * f / (k * t))


def sun(ws):
    ws.sunBlackbody()
    ws.sun.latitude = 20.0
    ws.sun.longitude = 35.0


def sun_direction():
    """Zenith and azimuth [deg] of the far sun seen from (lat, lon) = (0, 0)."""
    lat, lon = np.radians(20.0), np.radians(35.0)
    up, east, north = np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)
    return np.degrees(np.arccos(up)), np.degrees(np.arctan2(east, north))


def unit(zenith, azimuth):
    """A direction [north, east, up] from ARTS's zenith and azimuth (clockwise from north) [deg]."""
    z, a = np.radians(zenith), np.radians(azimuth)
    return np.stack(np.broadcast_arrays(np.sin(z) * np.cos(a), np.sin(z) * np.sin(a), np.cos(z)), axis=-1)


# %% 1. Clear sky against the exact layered solution

MIN_TAU = 1e-6
POLARIZATION = np.array([1.0, 0.3, -0.2, 0.1])  # the sky's [I, Q, U, V] / B(T_cmb)

ws = workspace([])
ws.legendre_polynomial_dimension = 1
ws.fourier_mode_dimension = 1
ws.vdisort_settings_agendaSetup(min_optical_depth=MIN_TAU)
ws.vdisort_settings_agendaExecute()

# An azimuth-independent Stokes vector is [I, Q, 0, 0] in the cosine and [0, 0, U, V] in the sine system of
# Fourier mode 0 (vdisort.cosine_mode, vdisort.sine_mode)
sky_boundary = np.asarray(ws.vdisort_settings.downward_boundary_condition)
cmb = sky_boundary[:, 0, 0, 0, 0].copy()
sky_boundary[:, 0, 0, :, :2] = cmb[:, None, None] * POLARIZATION[:2]
sky_boundary[:, 1, 0, :, 2:] = cmb[:, None, None] * POLARIZATION[2:]

ws.vdisort_spectral_rad_fieldCalc()
data = np.asarray(ws.vdisort_spectral_rad_field.data)
m = np.cos(np.radians(np.asarray(ws.vdisort_spectral_rad_field.zen_grid)))

ws.atm_pathFromPath()
ws.freq_grid_pathFromPath()
ws.spectral_propmat_pathFromPath()
temperature = np.array([p.temperature for p in ws.atm_path])
absorption = np.array([[pm[iv][0] for iv in range(len(FREQ))] for pm in ws.spectral_propmat_path])
assert (absorption[:2] == 0).all() and (absorption[2:] > 0).all(), "the top layer must be empty and no other"

for iv, f in enumerate(FREQ):
    dtau = np.maximum(-np.diff(Z) * 0.5 * (absorption[1:, iv] + absorption[:-1, iv]), MIN_TAU)
    tau = np.concatenate(([0.0], np.cumsum(dtau)))
    B = planck(f, temperature)
    c1 = np.diff(B) / dtau
    c0 = B[:-1] - c1 * tau[:-1]
    exact = np.zeros(data.shape[1:2] + (len(m), 4))
    for k, mk in enumerate(m):
        mu = abs(mk)
        e = np.exp(-dtau / mu)
        rise = -np.diff(B) * mu * np.expm1(-dtau / mu) / dtau  # (B_b - B_t) mu (1 - e) / dtau
        if mk < 0:  # downward from the sky
            I = cmb[iv] * POLARIZATION
            for l in range(len(dtau)):
                I = I * e[l]
                I[0] += B[l + 1] - B[l] * e[l] - rise[l]
                exact[l, k] = I
        else:  # upward from the blackbody surface; the field omits the top level
            I = np.array([planck(f, T_SURF), 0.0, 0.0, 0.0])
            exact[-1, k] = I
            for l in reversed(range(1, len(dtau))):
                I = I * e[l]
                I[0] += B[l] - B[l + 1] * e[l] + rise[l]
                exact[l - 1, k] = I
    largest = np.r_[np.abs(c0) + np.abs(c1) * (tau[1:] + 1.0), planck(f, T_SURF), cmb[iv]].max()
    tolerance = 16 * EPS * np.r_[largest, len(dtau) * np.abs(cmb[iv] * POLARIZATION[1:])]
    error = np.abs(data[iv, :, 0] - exact).max(axis=(0, 1))
    print(f"clear sky {f / 1e9:g} GHz: |VDISORT - exact| per Stokes {error}, tolerance {tolerance}")
    assert (error <= tolerance).all(), f"clear sky {f / 1e9:g} GHz: {error} > {tolerance}"
    transmitted = exact[-1, m < 0, 1] / (cmb[iv] * POLARIZATION[1])
    assert transmitted.min() < 1e-3 and transmitted.max() > 0.1, f"the streams must see the sky attenuated, {transmitted}"


# %% 2. Single scattering of the sun by an optically thin Henyey-Greenstein layer

EXTINCTION, ALBEDO, G = 1e-10, 0.9, 0.6
AZIMUTHS = [0.0, 37.0, 90.0, 140.0, 200.0, 270.0, 315.0]

ws = workspace(
    [A.HenyeyGreensteinScatterer(hg_ext, hg_ssa, G)],
    gas=False,
    fields={hg_ext: np.full_like(Z, EXTINCTION), hg_ssa: np.full_like(Z, ALBEDO)},
)
sun(ws)
ws.legendre_polynomial_dimension = LEGENDRE
ws.fourier_mode_dimension = LEGENDRE
ws.vdisort_settings_agendaSetup(
    layer_emission_setting="None",
    scattering_setting="ScatteringSpecies",
    space_setting="None",
    sun_setting="Sun",
    surf_setting="None",
)
data, m, azi = field(ws, AZIMUTHS)
s = ws.vdisort_settings

# The settings' sun is in the direction of the sun, to its parallax over the Earth's radius (below 0.0025 deg)
sza, saa = np.asarray(s.solar_zenith_angle), np.asarray(s.solar_azimuth_angle)
assert np.allclose(sza, sun_direction()[0], atol=2.5e-3) and np.allclose(saa, sun_direction()[1], atol=2.5e-3), (sza, saa)

legendre = (2 * np.arange(LEGENDRE) + 1) * G ** np.arange(LEGENDRE)
zenith = np.degrees(np.arccos(m))
for iv, f in enumerate(FREQ):
    tau = np.asarray(s.optical_thicknesses)[iv]
    assert np.allclose(np.asarray(s.single_scattering_albedo)[iv], ALBEDO, rtol=1e-14, atol=0)
    F0, mu0, total = np.asarray(s.solar_source)[iv, 0], np.cos(np.radians(sza[iv])), tau[-1]
    beam = unit(180.0 - sza[iv], saa[iv] + 180.0)
    cos_angle = unit(zenith[None, :], azi[:, None]) @ beam  # [azimuth, zenith]
    J = ALBEDO * F0 / (4 * np.pi) * np.polynomial.legendre.legval(cos_angle, legendre)
    mu = np.abs(m)
    t = tau[:, None, None]
    down = mu0 / (mu0 - mu) * np.exp(-t / mu0) * -np.expm1(-t * (1 / mu - 1 / mu0))
    up = mu0 / (mu0 + mu) * np.exp(-t / mu0) * -np.expm1(-(total - t) * (1 / mu0 + 1 / mu))
    first = J * np.where(m < 0, down, up)  # [layer, azimuth, zenith]

    grid = np.linspace(-1, 1, 20001)
    q = ALBEDO * total * np.abs(np.polynomial.legendre.legval(grid, legendre)).mean() / mu.min()
    scale = np.abs(first).max()
    error = np.abs(data[iv] - np.concatenate((first[..., None], np.zeros(first.shape + (3,))), axis=-1)).max() / scale
    print(f"single scattering {f / 1e9:g} GHz: |VDISORT - first order| / max = {error:.2e}, tolerance {q / (1 - q):.2e}")
    assert error <= q / (1 - q), f"single scattering {f / 1e9:g} GHz: {error} > {q / (1 - q)}"


# %% 3. Polarized scattering against RT3


def rayleigh():
    return A.GasScatterer(A.ConstantGasScattering(2e-31), A.RayleighGasScattering(0.0))


def henyey_greenstein():
    return A.HenyeyGreensteinScatterer(hg_ext, hg_ssa, 0.2)


def mie_drops():
    """Liquid drops as a Legendre series, converted from gridded Mie data."""
    x, _ = A.math.leggauss(4000)
    angles = np.degrees(np.arccos(np.asarray(x)[::-1]))
    t_grid = [220.0, 250.0, 280.0, 310.0]
    habit = A.ParticleHabit.liquid_sphere(t_grid, FREQ, [1.5e-3], A.IrregularZenithAngleGrid(angles))
    habit, reports = habit.to_tro_spectral_with_report(t_grid, FREQ, 64)
    for report in reports:
        assert np.abs(np.asarray(report.normalisation_error)).max() < 1e-7
    return A.ScatteringHabit(habit, A.MonodispersePSD(drops, 0.0, 400.0), 1.0, 3.0)


def against_rt3(name, ws, aziorder, azimuths, ground, thermal, minimum):
    ws.legendre_polynomial_dimension = LEGENDRE
    ws.fourier_mode_dimension = aziorder + 1
    data, m, azi = field(ws, azimuths)
    s = ws.vdisort_settings
    ws.atm_pathFromPath()
    ws.freq_grid_pathFromPath()
    ws.spectral_propmat_pathFromPath()

    settings = A.rt3.PathSettings(
        nstokes=4,
        nmu=NMU,
        quad=A.polradtran.QuadratureType.double_gauss,
        aziorder=aziorder,
        max_delta_tau=MAX_DELTA_TAU,
        normalisation_tolerance=1e-6,
    )
    h, k, c = A.constants.h, A.constants.k, A.constants.c
    for iv, f in enumerate(FREQ):
        # The sky temperature of the downward boundary, by inverting Planck's law
        sky = np.asarray(s.downward_boundary_condition)[iv, 0, 0, 0, 0]
        sky_temperature = h * f / k / np.log1p(2 * h * f**3 / (c**2 * sky)) if sky > 0 else 0.0
        problem = A.rt3.problem_from_path(
            ws.ray_path, ws.atm_path, ws.spectral_propmat_path, ws.freq_grid, iv, ws.scat_species,
            settings, ground, T_SURF, sky_temperature,
        )  # fmt: skip
        problem.thermal = thermal
        mu0 = 1.0
        if np.asarray(s.solar_source)[iv, 0] > 0:
            mu0 = np.cos(np.radians(np.asarray(s.solar_zenith_angle)[iv]))
            problem.direct_flux, problem.direct_mu = np.asarray(s.solar_source)[iv, 0] * mu0, mu0
        result = A.rt3.solve(problem)

        # VDISORT's phi is -azimuth and its beam propagates away from the sun; RT3's beam is at phi = 0
        phi = np.radians(np.asarray(s.solar_azimuth_angle)[iv] + 180.0 - azi)
        up, down = (np.asarray(A.rt3.azimuth_radiance(x, phi)) for x in (result.up, result.down))
        mu = np.asarray(result.mu)
        reference = np.empty_like(data[iv])
        for i, mi in enumerate(m):
            j = np.argmin(np.abs(mu - abs(mi)))
            assert abs(mu[j] - abs(mi)) < 1e-14, "the solvers' streams differ"
            reference[:, :, i] = (up if mi > 0 else down)[1:, :, j]

        scale = np.abs(reference[..., 0]).max()
        difference = np.abs(data[iv] - reference).max(axis=(0, 1, 2)) / scale
        polarization = np.abs(reference[..., 1:]).max(axis=(0, 1, 2)) / scale
        tolerance = 10 * MAX_DELTA_TAU / mu0
        print(f"{name} {f / 1e9:g} GHz: max |VDISORT - RT3| / max I per Stokes {difference}, tolerance {tolerance:.1e}")
        print(f"    max |Q|, |U|, |V| / max I = {polarization}")
        assert difference.max() <= tolerance, f"{name} {f / 1e9:g} GHz: {difference} > {tolerance}"
        assert (polarization >= minimum).all(), f"{name}: polarization {polarization} must exceed {minimum}"


ws = workspace([rayleigh(), mie_drops()])
ws.vdisort_settings_agendaSetup(scattering_setting="ScatteringSpecies", surf_setting="Thermal")
against_rt3("thermal", ws, 0, [0.0], A.polradtran.LambertianSurface(0.0), True, [1e-4, 0.0, 0.0])

ws = workspace([rayleigh(), henyey_greenstein(), mie_drops()])
sun(ws)
ws.vdisort_settings_agendaSetup(
    layer_emission_setting="None",
    scattering_setting="ScatteringSpecies",
    space_setting="None",
    sun_setting="Sun",
    surf_setting="Lambertian",
    surf_lambertian_value=[0.3] * len(FREQ),
)
against_rt3(
    "solar",
    ws,
    7,
    [0.0, 30.0, 75.0, 135.0, 180.0, 250.0],
    A.polradtran.LambertianSurface(0.3),
    False,
    [1e-2, 1e-2, 1e-5],
)
