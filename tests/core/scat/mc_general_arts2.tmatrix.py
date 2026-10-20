"""ARTS2 passive-MC I/Q acceptance cases using T-matrix-generated ARO ice data."""

from pathlib import Path

import numpy as np
import pyarts3 as pyarts
from scipy.interpolate import RegularGridInterpolator


A = pyarts.arts
DATA = Path(__file__).with_name("data") / "mc_general_arts2"


# Warren/Wiscombe/Gao REFICE, the material model used by historical arts_scat.
# Source: https://radiativetransfer.org/svn/rt/PyARTS/trunk/src/REFICE.f
# Evaluated with REAL inputs wavelength = c/f * 1e6 [micrometres], T [K],
# and COMPLEX (single precision) output, matching its original interface.
# Rows: frequency; columns: temperature. These are NOT fitted optical data.
REFRACTIVE_INDEX = np.array([
    [1.7800652980804443 + 0.00246472773142159j,
     1.7805060148239136 + 0.003070596605539322j],
    [1.7800800800323486 + 0.0024795890785753727j,
     1.7805144786834717 + 0.0030892009381204844j],
])
FREQUENCIES = [229e9, 231e9]
TEMPERATURES = [214.0, 225.0]
ANGLES = np.arange(0.0, 181.0, 10.0)
C = 299792458.0
# ARTS2 scattering data use Mishchenko's Stokes convention, whose U is the
# negative of ARTS's.  The legacy layout below holds them that way: D Z D with
# D = diag(1, 1, -1, 1) takes ARTS's phase matrices to ARTS2's.
D = np.diag([1.0, 1.0, -1.0, 1.0])


def generate(legacy_numerics=False, quadrature=10):
    # Match the historical integration in theta and relative azimuth, both
    # 0..180 degrees, with sin(theta) and a factor 2 for the other azimuth half.
    if legacy_numerics and quadrature == 10:
        # Original arts_math.py printed its quadrature to eight decimal places.
        x = np.array([.14887434, .43339539, .67940957, .86506337, .97390653])
        w = np.array([.29552422, .26926672, .21908636, .14945135, .06667134])
        x, w = np.r_[x, -x], np.r_[w, w]
    else:
        x, w = np.polynomial.legendre.leggauss(quadrature)
    theta = 90 * (1 + x)
    weights = 2 * (np.pi / 2)**2 * w[:, None] * \
        w[None, :] * np.sin(np.deg2rad(theta))[:, None]

    # The oblate symmetry axis stays vertical, alpha = beta = 0, and the
    # incident azimuth is 0.
    G = A.tmatrix.Geometry
    phase_geometries = [G(theta_incident=i, theta_scattered=s, phi_scattered=a)
                        for s in ANGLES for a in ANGLES for i in ANGLES]
    forward_geometries = [G(theta_incident=z, theta_scattered=z) for z in ANGLES]
    integral_geometries = [G(theta_incident=i, theta_scattered=s, phi_scattered=a)
                           for s in ANGLES for i in theta for a in theta]
    geometries = phase_geometries + forward_geometries + integral_geometries
    nphase = len(phase_geometries)
    nz = len(ANGLES)
    phase = np.empty((2, 2, nz, nz, nz, 1, 16))
    extinction = np.empty((2, 2, nz, 1, 3))
    absorption = np.empty((2, 2, nz, 1, 2))
    orders = []
    for fi, frequency in enumerate(FREQUENCIES):
        wavelength = C / frequency
        radius = 100e-6
        unit = 1.0
        if legacy_numerics:
            # arts_scat stored wavelength in float32 micrometres before calling
            # the double-precision solver, then converted amplitudes back to SI.
            wavelength = float(np.float32(wavelength * 1e6))
            radius, unit = 100., 1e-6
        for ti, _ in enumerate(TEMPERATURES):
            m = REFRACTIVE_INDEX[fi, ti]
            results = A.tmatrix.fixed_batch(radius, wavelength, A.tmatrix.Shape.spheroid(1.5), m.real, m.imag,
                                            geometries, accuracy=.001)
            orders.append(int(results[0].order))
            phase[fi, ti] = (np.array([D @ np.asarray(r.phase) @ D
                             for r in results[:nphase]]) * unit**2).reshape(nz, nz, nz, 1, 16)
            # Optical theorem: directional extinction, not the orientation-
            # averaged scalar cross section returned by the T-matrix solve.
            for zi, r in enumerate(results[nphase:nphase + nz]):
                s = np.asarray(r.amplitude) * unit
                lam = wavelength * unit
                k = A.Propmat([lam * (s[0, 0] + s[1, 1]).imag,
                               lam * (s[0, 0] - s[1, 1]).imag,
                               0., 0., 0., 0., lam * (s[1, 1] - s[0, 0]).real])
                extinction[fi, ti, zi, 0] = np.asarray(k)[[0, 1, 6]]
            # Integrate column 0 at fixed outgoing direction, as in arts_scat;
            # subtract scattering out of thermal equilibrium from K[:,0].
            z = np.array([D @ np.asarray(r.phase) @ D for r in results[nphase + nz:]]) * unit**2
            z = z.reshape(nz, quadrature, quadrature, 4, 4)
            scattering = np.einsum('ij,zijc->zc', weights, z[..., :, 0])
            for zi in range(nz):
                a = A.Stokvec([extinction[fi, ti, zi, 0, 0] - scattering[zi, 0],
                               extinction[fi, ti, zi, 0, 1] - scattering[zi, 1], 0., 0.])
                absorption[fi, ti, zi, 0] = np.asarray(a)[:2]
    ssd = A.SingleScatteringData()
    ssd.ptype = A.PType.PTYPE_AZIMUTH_RND
    ssd.description = (
        "ARTS3 T-matrix MC particle: equal-volume radius 100 um, oblate spheroid "
        "aspect ratio 1.5, alpha=beta=0, accuracy=0.001, original REFICE material. "
        f"Legacy numerics: {legacy_numerics}; absorption quadrature: {quadrature}. "
        "Full angular grid evaluated directly."
    )
    ssd.f_grid, ssd.T_grid = FREQUENCIES, TEMPERATURES
    ssd.za_grid, ssd.aa_grid = ANGLES, ANGLES
    ssd.pha_mat_data, ssd.ext_mat_data, ssd.abs_vec_data = phase, extinction, absorption
    return ssd, orders


def main():
    # Original ARTS2 grid literals, including their printed precision and the
    # distinct near-zero longitude node. No grid regeneration or rounding.
    p = np.array([
        101300.0, 90167.0235741, 80257.5729538, 71437.1813675, 63586.1600832,
        56597.9742864, 50377.7974505, 44841.2245837, 39913.1268917, 35526.6323136,
        31622.218098, 28146.9031067, 25053.5288841, 22300.1197385, 19849.3131507,
        18201.0, 18200.0, 17667.8527817, 17511.9573392, 16849.9258159,
        16212.922091, 15726.1371992, 15600.0, 13997.8181991, 12459.4432729,
        11090.1373673, 9871.31961931, 8786.45122233, 7820.81100195, 6961.29565628,
        6196.24195011, 5515.26845002, 4909.13465302, 4369.61559711, 3889.39025226,
        3461.94217734, 3081.47109493, 2742.8141842, 2441.37602376, 2173.06623384,
        1934.24397171, 1721.6685271, 1532.45534719, 1364.03689453, 1214.12780676,
        1080.69388524, 961.924491872, 856.207979617, 762.109823124, 678.353152889,
        603.801428709, 537.443017339, 478.377465095, 425.803279098, 379.007051376,
        337.353778244, 300.278243591, 267.277349148, 237.903287675, 211.757466418,
        188.485098386, 167.770388051, 149.332246144, 132.920475404, 118.312375511,
        105.309721143, 93.7360721512, 83.4343793426, 74.2648533966, 66.1030679856,
        58.838271366, 52.3718835274, 46.6161585058, 41.492993711, 36.9328701096,
        32.8739088828, 29.2610317593, 26.0452136273, 23.1828172865, 20.6350013108,
        18.3671929876, 16.3486191817, 14.5518887577, 12.952620895, 11.5291142506,
        10.262052482, 9.13424213291, 8.13037932606, 7.23684209634, 6.44150554692,
        5.73357731986, 5.10345114871, 4.54257650578, 4.04334257536, 3.59897497839,
        3.20344384717, 2.85138200283, 2.53801212505, 2.25908192607, 2.010806449,
        1.7898167077, 1.59311397114, 1.41802907199, 1.26218618719, 1.12347059917,
        1.0,
    ])
    lat = np.array([
        -25.0, -15.3469387755, -14.693877551, -14.0408163265, -13.387755102,
        -12.7346938776, -12.0816326531, -11.4285714286, -10.7755102041, -10.1224489796,
        -9.4693877551, -8.81632653061, -8.16326530612, -7.51020408163, -6.85714285714,
        -6.20408163265, -5.55102040816, -4.89795918367, -4.24489795918, -3.59183673469,
        -2.9387755102, -2.28571428571, -2.1, -2.0, -1.8,
        -1.63265306122, -1.08, -1.0, -0.979591836735, -0.36,
        -0.326530612245, 0.0, 0.326530612245, 0.36, 0.979591836735,
        1.0, 1.08, 1.63265306122, 1.8, 2.0,
        2.28571428571, 2.9387755102, 3.59183673469, 4.24489795918, 4.89795918367,
        5.55102040816, 6.20408163265, 6.85714285714, 7.51020408163, 8.16326530612,
        8.81632653061, 9.4693877551, 10.1224489796, 10.7755102041, 11.4285714286,
        12.0816326531, 12.7346938776, 13.387755102, 14.0408163265, 14.693877551,
        15.3469387755, 25.0,
    ])
    lon = np.array([
        -25.0, -15.3469387755, -14.693877551, -14.0408163265, -13.387755102,
        -12.7346938776, -12.0816326531, -11.4285714286, -10.7755102041, -10.1224489796,
        -9.4693877551, -8.81632653061, -8.16326530612, -7.51020408163, -6.85714285714,
        -6.20408163265, -5.55102040816, -4.89795918367, -4.24489795918, -3.59183673469,
        -2.9387755102, -2.28571428571, -2.1, -2.0, -1.8,
        -1.63265306122, -1.2, -1.0, -0.979591836735, -0.6,
        -0.326530612245, -2.22044604925e-16, 0.0, 0.326530612245, 0.6,
        0.979591836735, 1.0, 1.2, 1.63265306122, 1.8,
        2.0, 2.28571428571, 2.9387755102, 3.59183673469, 4.24489795918,
        4.89795918367, 5.55102040816, 6.20408163265, 6.85714285714, 7.51020408163,
        8.16326530612, 8.81632653061, 9.4693877551, 10.1224489796, 10.7755102041,
        11.4285714286, 12.0816326531, 12.7346938776, 13.387755102, 14.0408163265,
        14.693877551, 15.3469387755, 25.0,
    ])
    # Reconstruct AtmFieldsCalcExpand1D from source profiles, not its serialized
    # 3D outputs. ARTS2 interpolates these quantities linearly in log pressure.
    raw_p = A.GriddedField3.fromxml("planets/Earth/afgl/tropical/p.xml")
    raw_t = A.GriddedField3.fromxml("planets/Earth/afgl/tropical/t.xml")
    np.testing.assert_array_equal(raw_p.grids[0], raw_t.grids[0])
    log_p = np.log(np.asarray(raw_p.data).ravel()[::-1])
    z = np.interp(np.log(p), log_p, np.asarray(raw_p.grids[0])[::-1])
    t_profile = np.interp(np.log(p), log_p, np.asarray(raw_t.data).ravel()[::-1])
    shape = (len(p), len(lat), len(lon))
    t = np.broadcast_to(t_profile[:, None, None], shape)
    vmr = []
    for species in ("O2", "N2", "H2O"):
        raw = A.GriddedField3.fromxml(f"planets/Earth/afgl/tropical/{species}.xml")
        np.testing.assert_array_equal(raw_p.grids[0], raw.grids[0])
        profile = np.interp(
            np.log(p),
            log_p,
            np.asarray(raw.data).ravel()[::-1],
        )
        vmr.append(np.broadcast_to(profile[:, None, None], shape))

    # Original cloudbox indices: [13, 19, 23, 39, 23, 40], inclusive.
    # Interpolate the raw particle population in (log pressure, latitude, longitude).
    # Original pnd_field_raw.xml: a constant population inside a zero-padded box.
    pnd_pressure = [
        110000.0,
        21618.7922264,
        21617.7922264,
        20390.6782501,
        19233.2202635,
        18141.4643086,
        17111.6808705,
        17110.6808705,
        1e-5,
    ]
    pnd_latitude = [-90.0, -1.9, -1.8, -1.08, -0.36, 0.36, 1.08, 1.8, 1.9, 90.0]
    pnd_longitude = [
        -180.0,
        -1.9,
        -1.8,
        -1.2,
        -0.6,
        -2.22044604925e-16,
        0.6,
        1.2,
        1.8,
        1.9,
        180.0,
    ]
    raw_pnd = np.zeros((9, 10, 11))
    raw_pnd[2:7, 2:8, 2:9] = 34712.0922774  # particles / m³
    pnd_interpolate = RegularGridInterpolator(
        (np.log(pnd_pressure), pnd_latitude, pnd_longitude), raw_pnd
    )
    cloud_points = np.stack(
        np.meshgrid(np.log(p[13:20]), lat[23:40], lon[23:41], indexing="ij"), axis=-1
    )
    pnd = pnd_interpolate(cloud_points)

    ws = pyarts.Workspace()
    # Keep the original MC pressure-grid top, not the full AFGL 120-km extent.
    ws.atm_field.top_of_atmosphere = float(z[-1])
    grids = [z, lat, lon]
    names = ["altitude", "latitude", "longitude"]

    def field(data, field_grids=grids):
        return A.GeodeticField3(data=data, grids=field_grids, grid_names=names)

    ws.atm_field["t"] = field(t)
    ws.atm_field["p"] = field(np.broadcast_to(p[:, None, None], t.shape))
    for species, data in zip(("O2", "N2", "H2O"), vmr, strict=True):
        ws.atm_field[species] = field(data)

    number_density = A.ScatteringSpeciesProperty(
        "arts2_ice", A.ParticulateProperty.NumberDensity
    )
    pnd_full = np.zeros_like(t)
    pnd_full[13:20, 23:40, 23:41] = pnd
    ws.atm_field[number_density] = field(pnd_full)

    legacy, _ = generate()
    meta = A.ScatteringMetaData()
    # This is a one-particle population; size metadata only identifies its sole bin.
    meta.mass = 1.0
    meta.diameter_volume_equ = 1.0
    meta.diameter_max = 1.0
    habit = A.ParticleHabit.from_legacy_aro([legacy], [meta])
    psd = A.MonodispersePSD(number_density, 0.0, 400.0)
    ws.scat_species = [A.ScatteringHabit(habit, psd, 1.0, 1.0)]

    ws.abs_speciesSet(species=["O2-PWR98", "N2-SelfContStandardType", "H2O-PWR98"])
    # These predefined models have no fitted catalog payload; register their model
    # names locally so the test remains independent of tmp/ and ARTS_DATA_PATH.
    ws.abs_predef_dataReadSpeciesSplitCatalog(
        basename=str(DATA) + "/", name_missing=1, ignore_missing=1
    )
    ws.spectral_propmat_agendaAuto()
    ws.max_stepsize = 3e3
    ws.ray_path_observer_agendaSetGeometric(
        max_step_option="step", add_crossings=1, remove_non_atm=1
    )
    ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")
    ws.spectral_rad_surface_agendaSet(option="Blackbody")
    ws.surf_fieldPlanet(option="Earth")
    ws.surf_field[A.SurfaceKey("h")] = 500.0
    ws.surf_field[A.SurfaceKey("t")] = float(np.interp(500.0, z, t[:, 0, 0]))

    common = dict(
        frequency=230e9,
        sensor_pos=[95000.1, 7.61968838781, 0.0],
        sensor_los=[99.7841941981, 180.0],
        mc_seed=1729,
        mc_min_iter=128,
        mc_max_iter=128,
        mc_max_scatorder=30,
    )

    def check(reference):
        ws.MCGeneral(**common)
        rj = (
            np.asarray(ws.mc_spectral_rad)
            * A.constants.c**2
            / (2.0 * A.constants.k * common["frequency"] ** 2)
        )
        error = (
            np.asarray(ws.mc_error)
            * A.constants.c**2
            / (2.0 * A.constants.k * common["frequency"] ** 2)
        )
        accepted = np.abs(rj[:2] - reference) <= 4.0 * error[:2]
        assert np.all(accepted), (
            f"ARTS2 I/Q={np.asarray(reference)} K; "
            f"ARTS3 I/Q={rj[:2]} +/- {error[:2]} K; "
            f"ARTS3 U/V={rj[2:]} +/- {error[2:]} K"
        )

    ws.mc_antenna.set_pencil_beam()
    check((198.7, 7.9))

    ws.mc_antenna.set_gaussian_fwhm(0.1137, 0.239)
    check((198.6, 7.6))


if __name__ == "__main__":
    main()
