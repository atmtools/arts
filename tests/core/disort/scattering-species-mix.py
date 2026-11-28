"""DISORT's Legendre coefficients for a mix of scattering species.

DISORT takes the phase function of a layer as Legendre moments chi_l,
normalised to chi_0 = 1.  ARTS builds them in
disort_settingsLegendreCoefficientsFromPath from the bulk Legendre series of
all species at the layer's two levels, which only mixes the species right if
every species gives its series in the same normalisation (coefficients on the
orthonormal Y_l0, so that sqrt(4 pi) a_0 is its scattering coefficient,
K11 - a1, to the accuracy of its series).  The moments of the layer are then
the mean of each species' own normalised moments weighted by its scattering
coefficient, tabulated below from an independent evaluation of that mean.

The mix is liquid drops (a particle habit, converted to a Legendre series) and
Rayleigh-scattering gas, at two levels with different drop number densities.
"""

import numpy as np
import pyarts3 as pa

A = pa.arts
FREQ = 89e9
L = 16  # legendre_polynomial_dimension

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
    bulk = species.tro_spectral(point, A.AtmSpectralPoint(), [FREQ], L - 1)
    a = np.asarray(bulk.phase_matrix)[0, :, 0]
    sigma = np.asarray(bulk.extinction_matrix)[0, 0] - np.asarray(bulk.absorption_vector)[0, 0]
    return a, sigma, bulk.phase_matrix


# Every species' series is normalised alike: sqrt(4 pi) a_0 is its scattering coefficient
for point in levels:
    for species in (drop_species, gas_species):
        a, sigma, _ = series(species, point)
        assert abs(np.sqrt(4 * np.pi) * a[0] - sigma) < 1e-5 * sigma, (a[0], sigma)

# The scattering-weighted mean of the species' normalised moments over both levels
expected = np.array([
    1.0, 1.02706078840623e-01, 9.97614069928159e-02, 1.60120968410920e-02,
    4.38195762739369e-03, 4.12479214061834e-04, 2.39937783649114e-05, 9.73039350095837e-07,
    2.96845565583038e-08, 7.13218699449701e-10, 1.36522948827672e-11, 1.9374123442143e-13,
    0.0, 0.0, 0.0, 0.0,
])  # fmt: skip

# DISORT's moments of the mixed species from the workspace method
upper, lower = A.PropagationPathPoint(), A.PropagationPathPoint()
upper.pos = [1000.0, 0.0, 0.0]
lower.pos = [0.0, 0.0, 0.0]
ws = pa.Workspace()
ws.disort_settingsInit(
    freq_grid=[FREQ],
    ray_path=[upper, lower],
    quadrature_dimension=4,
    legendre_polynomial_dimension=L,
    fourier_mode_dimension=1,
)
ws.disort_settingsLegendreCoefficientsFromPath(
    spectral_phamat_spectral_path=[series(mix, point)[2] for point in levels]
)
got = np.asarray(ws.disort_settings.legendre_coefficients)[0, 0]
print("DISORT Legendre moments of drops + gas:", np.round(got[:6], 6))
np.testing.assert_allclose(got, expected, rtol=0, atol=1e-12)
