"""The VDISORT workspace methods against references derived outside VDISORT.

vdisort_settings_agendaSetup and vdisort_settings_agendaExecute build
VDISORT's input from ARTS data, and vdisort_spectral_rad_fieldCalc maps its
output onto the radiance field.  Three references check them:

1. Clear sky: ARTS's linear-in-tau transmission (rtepack.tran) and Planck
   function (physics.planck), layer by layer through the path's own
   propagation matrices and temperatures.  The sky is given a polarization
   to check Q, U and V; the blackbody surface keeps them zero upward.  The top
   layer has no absorption, so min_optical_depth sets its optical thickness.
   The solver writes the source as c0 + c1 tau in the global optical depth tau,
   so the tolerance of I is 16 eps times its largest term, |c0| + |c1| (tau + 1),
   or the largest boundary value.  Q, U and V are the sky's, transmitted layer
   by layer, so theirs is 16 eps per layer times the sky's value.

2. Single scattering of the sun by an optically thin Henyey-Greenstein layer:
   a table of the first-order radiance, per unit solar source, on two layers,
   all azimuths and four streams, and Q = U = V = 0 everywhere.  The sun's
   zenith and azimuth are tabulated too.  This pins the field's grids to
   propagation directions and the beam to the sun's position.  Higher orders
   are below q / (1 - q) of the largest first-order radiance, with
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

drops = A.ScatteringSpeciesProperty("drops", A.ParticulateProperty.NumberDensity)

ATMOSPHERE = {
    "t": 290.0 - 0.01 * Z,
    "p": 1e5 * np.exp(-Z / 8000.0),
    "H2O": np.where(GAS, 0.02 * np.exp(-Z / 2000.0), 0.0),
    "O2": np.where(GAS, 0.21, 0.0),
    drops: 30.0 * np.exp(-(((Z - 2500.0) / 1500.0) ** 2)),
}

# The spectral profiles of the Henyey-Greenstein scatterer, the same at all frequencies
HENYEY_GREENSTEIN = {
    "extinction": 2e-4 * np.exp(-(((Z - 3000.0) / 1000.0) ** 2)),
    "ssa": np.full_like(Z, 0.9),
    "g": np.full_like(Z, 0.2),
}


def spectral_profile(values):
    """A spectral field of a profile on Z, the same at all frequencies."""
    return A.SpectralGeodeticField4(
        data=np.asarray(values)[::-1][None, :, None, None],
        grid_names=["frequency", "altitude", "latitude", "longitude"],
        grids=[[FREQ[0]], Z[::-1], [0.0], [0.0]],
    )


def workspace(species, gas=True, spectral={}):
    """A workspace with the atmosphere (spectral replaces its spectral profiles) and a path down through its levels."""
    ws = pa.Workspace()
    ws.freq_grid = FREQ
    ws.abs_speciesSet(species=["H2O-PWR98", "O2-PWR98"] if gas else [])
    if gas:
        ws.ReadCatalogData()
    ws.spectral_propmat_agendaAuto()
    ws.surf_fieldPlanet(option="Earth")
    ws.surf_field[A.SurfaceKey("t")] = T_SURF
    profile = [A.AtmPoint() for _ in Z]
    for key, values in ATMOSPHERE.items():
        if isinstance(key, str) and key in ("H2O", "O2"):
            key = A.SpeciesEnum(key)
        for point, value in zip(profile, values):
            point[key] = value
    ws.atm_fieldFromProfile(atm_profile=profile[::-1], alt_grid=Z[::-1], altitude_extrapolation="Linear")
    for key, values in {**HENYEY_GREENSTEIN, **spectral}.items():
        ws.atm_spectral_field[key] = spectral_profile(values)
    ws.atm_scat_species = species
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


def sun(ws):
    ws.sunBlackbody()
    ws.sun.latitude = 20.0
    ws.sun.longitude = 35.0


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
absorption = np.array([[pm.K[iv][0] for iv in range(len(FREQ))] for pm in ws.spectral_propmat_path])
assert (absorption[:2] == 0).all() and (absorption[2:] > 0).all(), "the top layer must be empty and no other"

K = [pm.K for pm in ws.spectral_propmat_path]
dz = -np.diff(Z)


def evolve(I, iv, mu, start, end):
    """ARTS's linear-in-tau step through the layer between levels start and end, J = [B, 0, 0, 0]."""
    l = min(start, end)
    k1, k2 = K[l][iv], K[l + 1][iv]
    if l == 0:  # the empty top layer, min_optical_depth thick
        k1 = k2 = A.Propmat(MIN_TAU / dz[l])
    t = A.rtepack.tran(k1, k2, dz[l] / mu)
    return np.asarray(t()) @ (I - J[start]) + np.asarray(t.linsrc()) @ (J[start] - J[end]) + J[end]


for iv, f in enumerate(FREQ):
    B = np.asarray(A.physics.planck(f, temperature))
    J = np.zeros((len(Z), 4))
    J[:, 0] = B
    exact = np.zeros(data.shape[1:2] + (len(m), 4))
    for k, mk in enumerate(m):
        mu = abs(mk)
        if mk < 0:  # downward from the sky
            I = cmb[iv] * POLARIZATION
            for l in range(len(dz)):
                I = evolve(I, iv, mu, l, l + 1)
                exact[l, k] = I
        else:  # upward from the blackbody surface; the field omits the top level
            I = np.array([A.physics.planck(f, T_SURF), 0.0, 0.0, 0.0])
            exact[-1, k] = I
            for l in reversed(range(1, len(dz))):
                I = evolve(I, iv, mu, l + 1, l)
                exact[l - 1, k] = I
    # VDISORT's source c0 + c1 tau in its global optical depth
    tau = np.concatenate(([0.0], np.asarray(ws.vdisort_settings.optical_thicknesses)[iv]))
    c1 = np.diff(B) / np.diff(tau)
    c0 = B[:-1] - c1 * tau[:-1]
    largest = np.r_[np.abs(c0) + np.abs(c1) * (tau[1:] + 1.0), A.physics.planck(f, T_SURF), cmb[iv]].max()
    tolerance = 16 * EPS * np.r_[largest, len(dz) * np.abs(cmb[iv] * POLARIZATION[1:])]
    error = np.abs(data[iv, :, 0] - exact).max(axis=(0, 1))
    print(f"clear sky {f / 1e9:g} GHz: |VDISORT - exact| per Stokes {error}, tolerance {tolerance}")
    assert (error <= tolerance).all(), f"clear sky {f / 1e9:g} GHz: {error} > {tolerance}"
    transmitted = exact[-1, m < 0, 1] / (cmb[iv] * POLARIZATION[1])
    assert transmitted.min() < 1e-3 and transmitted.max() > 0.1, f"the streams must see the sky attenuated, {transmitted}"


# %% 2. Single scattering of the sun by an optically thin Henyey-Greenstein layer

EXTINCTION, ALBEDO, G = 1e-10, 0.9, 0.6
AZIMUTHS = [0.0, 37.0, 90.0, 140.0, 200.0, 270.0, 315.0]

ws = workspace(
    [A.HenyeyGreensteinScatterer("extinction", "ssa", "g")],
    gas=False,
    spectral={"extinction": np.full_like(Z, EXTINCTION), "ssa": np.full_like(Z, ALBEDO), "g": np.full_like(Z, G)},
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
data, m, _ = field(ws, AZIMUTHS)
s = ws.vdisort_settings

# The settings' sun is in the direction of the far sun at (lat, lon) = (20, 35) deg seen from (0, 0), zenith and
# azimuth (clockwise from north) in deg, to its parallax over the Earth's radius (below 0.0025 deg)
sza, saa = np.asarray(s.solar_zenith_angle), np.asarray(s.solar_azimuth_angle)
assert np.allclose(sza, 39.668454, atol=2.5e-3) and np.allclose(saa, 57.602310, atol=2.5e-3), (sza, saa)

# The first-order radiance / solar source on LAYERS, all AZIMUTHS and the STREAMS (zenith 40.3, 84.2, 95.8 and
# 139.7 deg), in units of 1e-6, and its largest value anywhere
LAYERS, STREAMS = [1, 5], [2, 6, 9, 13]
FIRST_ORDER = 1e-6 * np.array([
    [
        [0.010191, 0.107855, 0.041419, 0.013608],
        [0.009232, 0.084825, 0.032181, 0.010093],
        [0.009334, 0.087715, 0.034436, 0.010718],
        [0.011488, 0.137990, 0.056879, 0.019673],
        [0.016591, 0.347024, 0.166860, 0.088875],
        [0.016929, 0.373119, 0.182686, 0.103991],
        [0.013158, 0.186817, 0.077640, 0.029570],
    ],
    [
        [0.001699, 0.017976, 0.144965, 0.047629],
        [0.001539, 0.014137, 0.112632, 0.035326],
        [0.001556, 0.014619, 0.120526, 0.037514],
        [0.001915, 0.022998, 0.199077, 0.068857],
        [0.002765, 0.057837, 0.584009, 0.311064],
        [0.002821, 0.062187, 0.639400, 0.363970],
        [0.002193, 0.031136, 0.271739, 0.103496],
    ],
])  # fmt: skip
LARGEST = 3.187539e-6
MEAN_ABS_P = 1.0002  # mean of |P| over cos(angle) of the series truncated to LEGENDRE terms
assert np.isclose(np.degrees(np.arccos(m[STREAMS])), [40.29133, 84.16484, 95.83516, 139.70867], rtol=0, atol=1e-5).all()

for iv, f in enumerate(FREQ):
    tau = np.asarray(s.optical_thicknesses)[iv]
    assert np.allclose(np.asarray(s.single_scattering_albedo)[iv], ALBEDO, rtol=1e-14, atol=0)
    F0, total = np.asarray(s.solar_source)[iv, 0], tau[-1]

    q = ALBEDO * total * MEAN_ABS_P / np.abs(m).min()
    scale = F0 * LARGEST
    error = max(
        np.abs(data[iv][np.ix_(LAYERS, range(len(AZIMUTHS)), STREAMS)][..., 0] - F0 * FIRST_ORDER).max(),
        np.abs(data[iv][..., 1:]).max(),
    ) / scale
    print(f"single scattering {f / 1e9:g} GHz: |VDISORT - first order| / max = {error:.2e}, tolerance {q / (1 - q):.2e}")
    assert error <= q / (1 - q), f"single scattering {f / 1e9:g} GHz: {error} > {q / (1 - q)}"


# %% 3. Polarized scattering against RT3


def rayleigh():
    return A.GasScatterer(A.ConstantGasScattering(2e-31), A.RayleighGasScattering(0.0))


def henyey_greenstein():
    return A.HenyeyGreensteinScatterer("extinction", "ssa", "g")


def mie_drops():
    """Liquid drops as a Legendre series, converted from gridded Mie data."""
    t_grid = [220.0, 250.0, 280.0, 310.0]
    habit = A.ParticleHabit.liquid_sphere(t_grid, FREQ, [1.5e-3], A.GaussLegendreGrid(4000))
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
    for iv, f in enumerate(FREQ):
        # The sky temperature of the downward boundary
        sky = np.asarray(s.downward_boundary_condition)[iv, 0, 0, 0, 0]
        sky_temperature = A.physics.invplanck(sky, f) if sky > 0 else 0.0
        problem = A.rt3.problem_from_path(
            ws.ray_path, ws.atm_path, ws.atm_spectral_field, [pm.K for pm in ws.spectral_propmat_path], ws.freq_grid, iv,
            ws.atm_scat_species,
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
