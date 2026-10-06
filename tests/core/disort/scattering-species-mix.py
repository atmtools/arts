"""DISORT's Legendre coefficients for a mix of scattering species.

DISORT takes the phase function of a layer as Legendre moments chi_l,
normalised to chi_0 = 1.  ARTS builds them in
disort_settingsLegendreCoefficientsFromPath from the bulk Legendre series of
all species at the layer's two levels, which only mixes the species right if
every species gives its series in the same normalisation (coefficients on the
orthonormal Y_l0, so that sqrt(4 pi) a_0 is its scattering coefficient).  The
expected moments are then the mean of each species' own normalised moments
chi_s,l = a_s,l / (sqrt(2 l + 1) a_s,0), weighted by its scattering
coefficient sigma_s = sqrt(4 pi) a_s,0:

    chi_l = sum_{levels, s} sigma_s chi_s,l / sum_{levels, s} sigma_s.

sigma_s equals K11 - a1 of the species to the accuracy of its series.

The mix is liquid drops (a particle habit, converted to a Legendre series) and
Rayleigh-scattering gas, at two levels with different drop number densities.
"""

import numpy as np
import pyarts3 as pa

A = pa.arts
FREQ = 89e9
L = 16  # disort_legendre_polynomial_dimension

drops = A.ScatteringSpeciesProperty("drops", A.ParticulateProperty.NumberDensity)
habit = A.ParticleHabit.liquid_sphere([280.0], [FREQ], [1.5e-3], A.IrregularZenithAngleGrid(np.linspace(0, 180, 1801)))
habit, (report,) = habit.to_tro_spectral_with_report([280.0], [FREQ], 2 * L)
assert np.abs(np.asarray(report.normalisation_error)).max() < 1e-5

drop = A.ScatteringHabit(habit, A.MonodispersePSD(drops), 1.0, 3.0)
gas = A.GasScatterer(A.ConstantGasScattering(3e-29), A.RayleighGasScattering(0.0))
drop_species, gas_species, mix = A.ArrayOfScatteringSpecies(), A.ArrayOfScatteringSpecies(), A.ArrayOfScatteringSpecies()
drop_species.add(drop)
gas_species.add(gas)
mix.add(drop)
mix.add(gas)


def level(number_density):
    point = A.AtmPoint()
    point.temperature = 280.0
    point.pressure = 8e4
    point[drops] = number_density
    return point


levels = [level(20.0), level(200.0)]


def series(species, point):
    bulk = species.get_bulk_scattering_properties_tro_spectral(point, [FREQ], L - 1)
    a = np.asarray(bulk.phase_matrix)[0, :, 0, 0].real
    sigma = np.asarray(bulk.extinction_matrix)[0, 0] - np.asarray(bulk.absorption_vector)[0, 0]
    return a, sigma, bulk.phase_matrix


# Every species' series is normalised alike: sqrt(4 pi) a_0 is its scattering coefficient
expected, weight = np.zeros(L), 0.0
for point in levels:
    for species in (drop_species, gas_species):
        a, sigma, _ = series(species, point)
        phase_integral = np.sqrt(4 * np.pi) * a[0]
        assert abs(phase_integral - sigma) < 1e-5 * sigma, (phase_integral, sigma)
        expected += phase_integral * a / (np.sqrt(2 * np.arange(L) + 1) * a[0])
        weight += phase_integral
expected /= weight
assert expected[1] > 0.1, "the drops must make the mix forward scattering"

# DISORT's moments of the mixed species from the workspace method
upper, lower = A.PropagationPathPoint(), A.PropagationPathPoint()
upper.pos = [1000.0, 0.0, 0.0]
lower.pos = [0.0, 0.0, 0.0]
ws = pa.Workspace()
ws.disort_settingsInit(
    freq_grid=[FREQ],
    ray_path=[upper, lower],
    disort_quadrature_dimension=4,
    disort_legendre_polynomial_dimension=L,
    disort_fourier_mode_dimension=1,
)
ws.disort_settingsLegendreCoefficientsFromPath(
    spectral_phamat_spectral_path=[series(mix, point)[2] for point in levels]
)
got = np.asarray(ws.disort_settings.legendre_coefficients)[0, 0]
print("DISORT Legendre moments of drops + gas:", np.round(got[:6], 6))
np.testing.assert_allclose(got, expected, rtol=0, atol=1e-12)
