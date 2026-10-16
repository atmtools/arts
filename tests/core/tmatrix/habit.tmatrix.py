"""T-matrix particles through the native habit/PSD/scattering-species path."""

import numpy as np
import pyarts3 as pa

A = pa.arts
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
        aspect_ratio=1.5,
        angles=NANGLE,
    )
    args.update(kw)
    return A.ParticleHabit.tmatrix(**args)


habit = make()
ssd = habit[0]
phase = np.array(ssd.phase_matrix)
extinction = np.array(ssd.extinction_matrix)
absorption = np.array(ssd.absorption_vector)
assert phase.shape == (2, 2, NANGLE, 6)
np.testing.assert_allclose(ssd.properties.mass, RHO * np.pi / 6 * D**3)
np.testing.assert_allclose(ssd.properties.d_veq, D)
np.testing.assert_allclose(ssd.properties.d_max, D * 1.5 ** (1 / 3))
for it in range(2):
    for jf in range(2):
        direct = A.tmatrix.random(
            D / 2, C / F[jf], 1.5, M[it, jf].real, M[it, jf].imag, angles=NANGLE
        )
        full = np.asarray(direct.phase)
        compact = full[:, [0, 0, 1, 2, 2, 3], [0, 1, 1, 2, 3, 3]]
        np.testing.assert_allclose(
            phase[it, jf], compact * direct.scattering / (4 * np.pi), rtol=1e-12
        )
        np.testing.assert_allclose(extinction[it, jf, 0], direct.extinction, rtol=1e-12)
        np.testing.assert_allclose(
            absorption[it, jf, 0], direct.extinction - direct.scattering, rtol=1e-12
        )
np.testing.assert_allclose(ssd.backscatter_matrix, phase[:, :, -1, :])
np.testing.assert_allclose(ssd.forwardscatter_matrix, phase[:, :, 0, :])

# Material density changes mass, not the optical response for a supplied index.
heavy = make(density=2 * RHO)
np.testing.assert_allclose(heavy[0].properties.mass, 2 * ssd.properties.mass)
np.testing.assert_array_equal(heavy[0].phase_matrix, phase)
multi = make(diameters=[D, 2 * D], shape=-2)
assert len(multi) == 2
length = D * (2 / (3 * 1.5**2)) ** (1 / 3)
np.testing.assert_allclose(multi[0].properties.d_max, np.hypot(length, 1.5 * length))
np.testing.assert_allclose(multi[1].properties.mass, 8 * multi[0].properties.mass)
assert np.max(np.abs(np.asarray(multi[0].phase_matrix) - phase)) > 1e-3 * np.max(
    np.abs(phase)
)

# Existing species machinery interpolates temperature/frequency and applies N.
number = A.ScatteringSpeciesProperty("tmatrix_ice", A.ParticulateProperty.NumberDensity)
psd = A.MonodispersePSD(number, 0.0, 400.0)
species = A.ArrayOfScatteringSpecies()
species.add(A.ScatteringHabit(habit, psd, RHO * np.pi / 6, 3.0))
point = A.AtmPoint()
point.temperature = 270.0
za = A.IrregularZenithAngleGrid(np.linspace(0, 180, NANGLE))
for concentration in [1e6, 2e6, 0, 1e6]:
    point[number] = float(concentration)
    bulk = species.get_bulk_scattering_properties_tro_gridded(point, [110e9], za)
    np.testing.assert_allclose(
        np.asarray(bulk.phase_matrix)[0, 0],
        concentration * phase.mean(axis=(0, 1)),
        rtol=1e-12,
        atol=1e-30,
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

# Independent small-sphere Rayleigh limit checks SI units and normalization.
r = 1e-6
m = 1.5 + 0.01j
ray = make(
    t_grid=[270.0],
    f_grid=[100e9],
    diameters=[2 * r],
    refractive_index=[[m]],
    aspect_ratio=1.0,
    angles=181,
)[0]
k = 2 * np.pi * 100e9 / C
contrast = (m * m - 1) / (m * m + 2)
csca = 8 * np.pi / 3 * k**4 * r**6 * abs(contrast) ** 2
cabs = 4 * np.pi * k * r**3 * contrast.imag
np.testing.assert_allclose(np.asarray(ray.absorption_vector)[0, 0, 0], cabs, rtol=1e-4)
mu = np.cos(np.linspace(0, np.pi, 181))
expected = csca / (4 * np.pi) * 0.75 * (1 + mu * mu)
np.testing.assert_allclose(
    np.asarray(ray.phase_matrix)[0, 0, :, 0], expected, rtol=1e-4
)

for kw in [
    dict(density=0),
    dict(diameters=[]),
    dict(t_grid=[280, 260]),
    dict(refractive_index=[[m]]),
    dict(refractive_index=np.full((2, 2), np.nan)),
    dict(angles=1),
    dict(shape=0),
]:
    try:
        make(**kw)
    except ValueError:
        pass
    else:
        raise AssertionError(f"Accepted invalid input {kw}")
