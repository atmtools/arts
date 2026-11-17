"""T-matrix particles through the native habit/PSD/scattering-species path."""

import numpy as np
import pyarts3 as pa

A = pa.arts
# Random orientation computes Clebsch-Gordan coefficients with the Wigner tables
pa.Workspace().WignerInit()
C = 299792458.0
T = [260.0, 280.0]
F = [100e9, 120e9]
M = np.array([[1.5 + 0.01j, 1.6 + 0.02j], [1.55 + 0.015j, 1.65 + 0.025j]])
D = 20e-6
RHO = 917.0
NANGLE = 19


def make(**kw):
    args = dict(
        t_grid=T,
        f_grid=F,
        diameters=[D],
        refractive_index=M,
        density=RHO,
        aspect_ratio=[1.5],
    )
    args.update(kw)
    return A.ParticleHabit.tmatrix(**args)


habit = make()
ssd = habit[0]
# The particle holds the exact series, to the degree of the longest
phase = np.array(ssd.phase_matrix)
extinction = np.array(ssd.extinction_matrix)
absorption = np.array(ssd.absorption_vector)
assert phase.shape[:2] == (2, 2) and phase.shape[3] == 6
ANGLES = np.linspace(0, 180, NANGLE)
# The series evaluated by ARTS at the angles
gridded = np.asarray(ssd.phase_matrix.to_gridded(A.IrregularZenithAngleGrid(ANGLES)))
np.testing.assert_allclose(ssd.properties.mass, RHO * np.pi / 6 * D**3)
np.testing.assert_allclose(ssd.properties.d_veq, D)
# The longer axis of the spheroid
np.testing.assert_allclose(ssd.properties.d_max, 2.289428485106664e-05)
for it in range(2):
    for jf in range(2):
        direct = A.tmatrix.random(
            D / 2, C / F[jf], A.tmatrix.Shape.spheroid(1.5), M[it, jf].real, M[it, jf].imag, angles=NANGLE
        )
        series = np.asarray(direct.legendre) * direct.scattering / (4 * np.pi)
        n = series.shape[0]
        np.testing.assert_allclose(phase[it, jf, :n], series, rtol=1e-12, atol=0)
        assert not phase[it, jf, n:].any()
        # The series is the scattering matrix the solver evaluates directly
        np.testing.assert_allclose(
            gridded[it, jf],
            np.asarray(direct.phase) * direct.scattering / (4 * np.pi),
            rtol=0,
            atol=1e-12 * np.abs(series).max(),
        )
        np.testing.assert_allclose(extinction[it, jf, 0], direct.extinction, rtol=1e-12)
        np.testing.assert_allclose(
            absorption[it, jf, 0], direct.extinction - direct.scattering, rtol=1e-12
        )
# The endpoint matrices are the series at 180 and 0 degrees, where |P_l| = 1.  The endpoints and the evaluation of
# the series may sum its L + 1 terms in their own order: the difference is the rounding of both sums, at most
# (L + 1) ulp each of the sum of the terms' magnitudes.  F12 and F34 vanish there, so they are that rounding alone.
degrees = np.arange(phase.shape[-2])
terms = np.einsum("...lk,l->...k", np.abs(phase), np.sqrt((2 * degrees + 1) / (4 * np.pi)))
rounding = 2 * (degrees.size + 1) * np.finfo(float).eps * terms
endpoints = np.asarray(ssd.phase_matrix.to_gridded(A.IrregularZenithAngleGrid([0.0, 180.0])))
for endpoint, k in [(ssd.backscatter_matrix, 1), (ssd.forwardscatter_matrix, 0)]:
    difference = np.abs(np.asarray(endpoint) - endpoints[:, :, k])
    assert np.all(difference <= rounding), (k, difference / rounding)

# Material density changes mass, not the optical response for a supplied index.
heavy = make(density=2 * RHO)
np.testing.assert_allclose(heavy[0].properties.mass, 2 * ssd.properties.mass)
np.testing.assert_array_equal(heavy[0].phase_matrix, phase)
multi = make(diameters=[D, 2 * D], shape="cylinder")
assert len(multi) == 2
# The diagonal of the cylinder of diameter 1.5 times its length
np.testing.assert_allclose(multi[0].properties.d_max, 2.403700850309326e-05)
np.testing.assert_allclose(multi[1].properties.mass, 8 * multi[0].properties.mass)
cylinder = np.asarray(multi[0].phase_matrix.to_gridded(A.IrregularZenithAngleGrid(ANGLES)))
assert np.max(np.abs(cylinder - gridded)) > 1e-3 * np.max(np.abs(gridded))

# A Chebyshev particle r_0 (1 + eps cos(n theta)) of odd degree, without mirror symmetry
n, eps = 3, 0.1
cheb = make(diameters=[D, 2 * D], shape="chebyshev", degree=n, aspect_ratio=[eps])
assert cheb[0].properties.name == "T-matrix Chebyshev particle"
np.testing.assert_allclose(cheb[1].properties.d_max, 2 * cheb[0].properties.d_max)
# Its maximum dimension against the largest distance of two surface points among 4001 polar angles, at the
# volume-equivalent diameter D
grid = 2.061498031474102e-05
assert grid <= cheb[0].properties.d_max <= grid * (1 + 1e-6)

# Per-diameter shapes: each particle is the one of its own entry
pair = make(diameters=[D, 2 * D], aspect_ratio=[1.5, 0.5])
for particle, single in [(pair[0], make()[0]), (pair[1], make(diameters=[2 * D], aspect_ratio=[0.5])[0])]:
    for field in ["phase_matrix", "extinction_matrix", "absorption_vector"]:
        np.testing.assert_array_equal(np.array(getattr(particle, field)), np.array(getattr(single, field)))
    assert particle.properties.d_max == single.properties.d_max

# A generalized Chebyshev particle r_0 (1 + c_0 + sum_n c_n cos(n theta)), whose coefficients are test data
DROP = [-0.0481, 0.0359, -0.1263, 0.0244, 0.0091, -0.0099, 0.0015, 0.0025, -0.0016, -0.0002, 0.0010]
drop = make(shape="generalized_chebyshev", coefficients=[DROP])
assert drop[0].properties.name == "T-matrix generalized Chebyshev particle"
# Its maximum dimension against the largest distance of two surface points among 4001 polar angles, at the
# volume-equivalent diameter D
grid = 2.17370105873609e-05
assert grid <= drop[0].properties.d_max <= grid * (1 + 1e-6)
direct = A.tmatrix.random(
    D / 2, C / F[0], A.tmatrix.Shape.generalized_chebyshev(DROP), M[0, 0].real, M[0, 0].imag
)
np.testing.assert_allclose(np.array(drop[0].extinction_matrix)[0, 0, 0], direct.extinction, rtol=1e-15)

# Existing species machinery interpolates temperature/frequency and applies N.
number = A.ScatteringSpeciesProperty("tmatrix_ice", A.ParticulateProperty.NumberDensity)
psd = A.MonodispersePSD(number, 0.0, 400.0)
species = A.ArrayOfScatteringSpecies()
species.add(A.ScatteringHabit(habit, psd, RHO * np.pi / 6, 3.0))
point = A.AtmPoint()
point.temperature = 270.0
za = A.IrregularZenithAngleGrid(ANGLES)
for concentration in [1e6, 2e6, 0, 1e6]:
    point[number] = float(concentration)
    bulk = species.get_bulk_scattering_properties_tro_gridded(point, [110e9], za)
    # The series of the habit's grid points, which surround the point equally; its rounding is relative to the
    # largest element
    expected = concentration * gridded.mean(axis=(0, 1))
    np.testing.assert_allclose(
        np.asarray(bulk.phase_matrix)[0, 0],
        expected,
        rtol=1e-12,
        atol=1e-13 * np.abs(expected).max(),
    )
    np.testing.assert_allclose(
        np.asarray(bulk.extinction_matrix)[0, 0],
        concentration * extinction.mean(axis=(0, 1)),
        rtol=1e-12,
        atol=1e-30,
    )
    np.testing.assert_allclose(
        np.asarray(bulk.absorption_vector)[0, 0],
        concentration * absorption.mean(axis=(0, 1)),
        rtol=1e-12,
        atol=1e-30,
    )
np.testing.assert_array_equal(habit[0].phase_matrix, phase)

# Solvers requesting laboratory-frame ARO data can consume a TRO habit too.
point[number] = 1e6
lab = species.get_bulk_scattering_properties_aro_gridded(
    point, [110e9], [0.0, 90.0, 180.0], [0.0, 90.0, 180.0], za
)
lab_phase = np.array(lab.phase_matrix)
assert np.isfinite(lab_phase).all()
point[number] = 2e6
lab_twice = species.get_bulk_scattering_properties_aro_gridded(
    point, [110e9], [0.0, 90.0, 180.0], [0.0, 90.0, 180.0], za
)
np.testing.assert_allclose(
    lab_twice.phase_matrix, 2 * lab_phase, rtol=1e-12, atol=1e-30
)
np.testing.assert_array_equal(habit[0].phase_matrix, phase)

# The small-sphere Rayleigh limit checks SI units and normalization: the absorption cross section
# 4 pi k r^3 Im((m^2 - 1) / (m^2 + 2)) and F11 = csca / (4 pi) 3 / 4 (1 + cos^2 Theta) at 0, 45, 90, 135 and 180 deg
r = 1e-6
m = 1.5 + 0.01j
ray = make(
    t_grid=[270.0],
    f_grid=[100e9],
    diameters=[2 * r],
    refractive_index=[[m]],
    aspect_ratio=[1.0],
)[0]
np.testing.assert_allclose(np.asarray(ray.absorption_vector)[0, 0, 0], 1.3122980284273886e-16, rtol=1e-4)
rayleigh = np.asarray(ray.phase_matrix.to_gridded(A.IrregularZenithAngleGrid([0.0, 45.0, 90.0, 135.0, 180.0])))
np.testing.assert_allclose(
    rayleigh[0, 0, :, 0],
    [1.6697761628094516e-24, 1.2523321221070887e-24, 8.348880814047258e-25, 1.2523321221070887e-24, 1.6697761628094516e-24],
    rtol=1e-4,
)

for kw in [
    dict(density=0),
    dict(diameters=[]),
    dict(t_grid=[280, 260]),
    dict(refractive_index=[[m]]),
    dict(refractive_index=np.full((2, 2), np.nan)),
    dict(shape="chebyshev"),
    dict(shape="chebyshev", degree=2, aspect_ratio=[1.5]),
    dict(degree=2),
    dict(aspect_ratio=[1.5, 2.0]),
    dict(aspect_ratio=[]),
    dict(coefficients=[[0.0, 0.1]]),
    dict(shape="generalized_chebyshev"),
    dict(shape="generalized_chebyshev", coefficients=[[0.0, 0.1], [0.0, 0.2]]),
    dict(shape="generalized_chebyshev", coefficients=[[0.0, 1.0]]),
]:
    try:
        make(**kw)
    except ValueError:
        pass
    else:
        raise AssertionError(f"Accepted invalid input {kw}")

# Azimuthally random orientation: every particle is tmatrix.azimuthally_random on its grids, in SI units
ZA_INC, DAA, ZA_SCAT = [0.0, 40.0, 90.0, 140.0, 180.0], [-180.0, -60.0, 0.0, 60.0, 180.0], [0.0, 40.0, 90.0, 140.0, 180.0]
TILTS, WEIGHTS = [10.0, 30.0], [2.0, 1.0]


def make_aro(**kw):
    args = dict(
        t_grid=T,
        f_grid=F,
        diameters=[D, 2 * D],
        refractive_index=M,
        density=RHO,
        aspect_ratio=[1.5, 0.5],
        za_inc_grid=ZA_INC,
        delta_aa_grid=DAA,
        za_scat_grid=A.IrregularZenithAngleGrid(ZA_SCAT),
        tilt_angles=TILTS,
        tilt_weights=WEIGHTS,
    )
    args.update(kw)
    return A.ParticleHabit.tmatrix_aro(**args)


aro = make_aro()
for particle, diameter, aspect, d_max in zip(aro, [D, 2 * D], [1.5, 0.5], [2.289428485106664e-05, 6.349604207872798e-05]):
    np.testing.assert_allclose(particle.properties.mass, RHO * np.pi / 6 * diameter**3)
    # The longer axis of the spheroid of horizontal to rotational axis ratio aspect
    np.testing.assert_allclose(particle.properties.d_max, d_max)
    phase = np.array(particle.phase_matrix)
    assert phase.shape == (2, 2, 5, 5, 5, 4, 4)
    for it in range(2):
        for jf in range(2):
            direct = A.tmatrix.azimuthally_random(
                diameter / 2,
                C / F[jf],
                A.tmatrix.Shape.spheroid(aspect),
                M[it, jf].real,
                M[it, jf].imag,
                ZA_INC,
                DAA,
                ZA_SCAT,
                TILTS,
                WEIGHTS,
            )
            np.testing.assert_array_equal(phase[it, jf], np.array(direct.phase))
            np.testing.assert_array_equal(np.array(particle.extinction_matrix)[it, jf], np.array(direct.extinction))
            np.testing.assert_array_equal(np.array(particle.absorption_vector)[it, jf], np.array(direct.absorption))
    # The grids hold the backward and forward directions, so the extracted matrices are grid values
    for k, za in enumerate(ZA_INC):
        np.testing.assert_array_equal(
            np.array(particle.backscatter_matrix)[:, :, k], phase[:, :, k, -1, ZA_SCAT.index(180.0 - za)]
        )
        np.testing.assert_array_equal(np.array(particle.forwardscatter_matrix)[:, :, k], phase[:, :, k, 2, k])

# A sphere's ARO habit is its TRO habit in the laboratory frame, which the Legendre series give exactly
sphere_args = dict(diameters=[D], aspect_ratio=[1.0])
bulk = []
for sphere in [make(**sphere_args), make_aro(**sphere_args)]:
    species = A.ArrayOfScatteringSpecies()
    species.add(A.ScatteringHabit(sphere, psd, RHO * np.pi / 6, 3.0))
    point[number] = 1e6
    bulk.append(
        species.get_bulk_scattering_properties_aro_gridded(
            point, [110e9], ZA_INC, DAA, A.IrregularZenithAngleGrid(ZA_SCAT)
        )
    )
tro_phase, aro_phase = np.array(bulk[0].phase_matrix), np.array(bulk[1].phase_matrix)
assert np.max(np.abs(aro_phase - tro_phase)) <= 1e-3 * np.max(np.abs(tro_phase))
np.testing.assert_allclose(
    np.array(bulk[1].extinction_matrix)[..., 0], np.array(bulk[0].extinction_matrix)[..., 0], rtol=1e-3
)

for kw in [
    dict(tilt_weights=[1.0]),
    dict(tilt_angles=[200.0, 10.0]),
    dict(za_inc_grid=[-1.0]),
    dict(za_inc_grid=[90.0, 40.0]),
    dict(delta_aa_grid=[0.0, 0.0]),
    dict(aspect_ratio=[1.5, 0.5, 1.0]),
    dict(density=0),
]:
    try:
        make_aro(**kw)
    except ValueError:
        pass
    else:
        raise AssertionError(f"Accepted invalid ARO input {kw}")
