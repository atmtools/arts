"""The spectral forms of scattering data through pyarts3: Legendre series and azimuthal Fourier modes.

- Gridded TRO data convert to a Legendre series with a report on how well it
  represents them; the series evaluates exactly on any grid.
- A habit of gridded data refuses to give a Legendre series, and names the
  conversion; converted, it gives one of the degree asked for.
- The azimuthal Fourier modes of the laboratory-frame phase matrix
  (get_bulk_scattering_properties_aro_fourier) evaluate to the gridded
  laboratory-frame phase matrix at any azimuth, exactly for Rayleigh
  scattering, whose modes above m = 2 vanish.
- A user's function of the Legendre series gives the Fourier modes of its
  series to the degree set on the species, and refuses without one.
- Beyond its temperature grid a habit extrapolates linearly.
"""

import numpy as np
import pyarts3 as pa

A = pa.arts

FREQ = 89e9
T = 280.0

# Gridded Mie data of 0.8 mm drops, a Legendre series of degree 32 and its report
angles = np.linspace(0.0, 180.0, 1801)
grid = A.IrregularZenithAngleGrid(angles)
habit = A.ParticleHabit.liquid_sphere([T], [FREQ], [0.8e-3], grid)
ssd = habit[0]
spectral, report = ssd.to_spectral_with_report(32)
assert spectral.phase_matrix.get_degree() == 32
assert np.abs(np.asarray(report.normalisation_error)).max() < 1e-5, np.asarray(report.normalisation_error)
assert np.asarray(report.tail).max() < 1e-10, "the series of these drops has converged at degree 32"
assert np.asarray(report.reconstruction_error).max() < 1e-5, np.asarray(report.reconstruction_error)
assert np.asarray(report.min_f11).min() > 0.0
gridded = np.asarray(ssd.phase_matrix)
series = np.asarray(spectral.phase_matrix.to_gridded(grid))
assert np.abs(series - gridded).max() <= np.asarray(report.reconstruction_error).max() * np.abs(gridded[..., 0]).max() * (
    1 + 1e-9
)

# A habit of gridded data has no Legendre series; converted, it has
drops = A.ScatteringSpeciesProperty("drops", A.ParticulateProperty.NumberDensity)
point = A.AtmPoint()
point.temperature = T
point.pressure = 1e5
point[drops] = 100.0
gridded_habit = A.ScatteringHabit(habit, A.MonodispersePSD(drops), 1.0, 3.0)
try:
    gridded_habit.get_bulk_scattering_properties_tro_spectral(point, [FREQ], 8)
    raise AssertionError("gridded habit data must not give a Legendre series")
except RuntimeError:
    pass
converted, reports = habit.to_tro_spectral_with_report([T], [FREQ], 32)
bulk = A.ScatteringHabit(converted, A.MonodispersePSD(drops), 1.0, 3.0).get_bulk_scattering_properties_tro_spectral(
    point, [FREQ], 8
)
assert np.asarray(bulk.phase_matrix).shape[:2] == (1, 9)
try:
    A.ScatteringHabit(converted, A.MonodispersePSD(drops), 1.0, 3.0).get_bulk_scattering_properties_tro_spectral(
        point, [FREQ], 40
    )
    raise AssertionError("a series of degree 32 must not give degree 40")
except RuntimeError:
    pass

# Fourier modes of Rayleigh scattering against the gridded laboratory frame
species = A.ArrayOfScatteringSpecies()
species.add(A.GasScatterer(A.ConstantGasScattering(1e-30), A.RayleighGasScattering(0.0)))
za_inc, za_scat, delta = [0.0, 30.0, 150.0], [20.0, 30.0, 160.0, 180.0], [0.0, 45.0, 180.0, 200.0, 359.0]
modes = species.get_bulk_scattering_properties_aro_fourier(point, [FREQ], za_inc, za_scat, 4)
lab = species.get_bulk_scattering_properties_aro_gridded(point, [FREQ], za_inc, delta, A.IrregularZenithAngleGrid(za_scat))
from_modes = np.asarray(modes.phase_matrix.to_gridded(delta))
reference = np.asarray(lab.phase_matrix)
assert np.abs(from_modes - reference).max() < 1e-12 * np.abs(reference).max(), np.abs(from_modes - reference).max()
assert np.abs(np.asarray(modes.phase_matrix)[:, :, :, :, 3:]).max() < 1e-13 * np.abs(reference).max()
try:
    modes.phase_matrix.to_fourier(5)
    raise AssertionError("modes to m = 4 must not give m = 5")
except RuntimeError:
    pass

# A user's function must give exactly the Legendre series asked for
rayleigh = A.GasScatterer(A.ConstantGasScattering(1e-30), A.RayleighGasScattering(0.0))
good, bad = A.ArrayOfScatteringSpecies(), A.ArrayOfScatteringSpecies()
good.add(A.ScatteringGeneralSpectralTRO(lambda atm, f, l: rayleigh.get_bulk_scattering_properties_tro_spectral(atm, f, l)))
bad.add(A.ScatteringGeneralSpectralTRO(lambda atm, f, l: rayleigh.get_bulk_scattering_properties_tro_spectral(atm, f, l + 1)))
user = good.get_bulk_scattering_properties_tro_spectral(point, [FREQ], 4)
builtin = rayleigh.get_bulk_scattering_properties_tro_spectral(point, [FREQ], 4)
assert np.array_equal(np.asarray(user.phase_matrix), np.asarray(builtin.phase_matrix))
try:
    bad.get_bulk_scattering_properties_tro_spectral(point, [FREQ], 4)
    raise AssertionError("a user function giving the wrong degree must fail")
except RuntimeError:
    pass

# The Fourier modes of a user's function are those of its series to the degree set on the species; the Rayleigh
# phase matrix is a polynomial of degree 2 in cos(Theta), so its series to degree 2 is the matrix itself
try:
    good.get_bulk_scattering_properties_aro_fourier(point, [FREQ], za_inc, za_scat, 4)
    raise AssertionError("the Fourier modes of a user function need the degree of its series")
except RuntimeError:
    pass
user_species = A.ArrayOfScatteringSpecies()
user_species.add(A.ScatteringGeneralSpectralTRO(lambda atm, f, l: rayleigh.get_bulk_scattering_properties_tro_spectral(atm, f, l), 2))
user_modes = user_species.get_bulk_scattering_properties_aro_fourier(point, [FREQ], za_inc, za_scat, 4)
reference = np.asarray(modes.phase_matrix)
difference = np.abs(np.asarray(user_modes.phase_matrix) - reference).max()
assert difference < 1e-12 * np.abs(reference).max(), difference
assert np.array_equal(np.asarray(user_modes.extinction_matrix), np.asarray(modes.extinction_matrix))
assert np.array_equal(np.asarray(user_modes.absorption_vector), np.asarray(modes.absorption_vector))


# Beyond its temperature grid a habit extrapolates linearly.  The regrids once dropped the negative weight there,
# which gave 1.6 times the end value at these temperatures.
two_temperatures = A.ParticleHabit.liquid_sphere(
    [250.0, 300.0], [FREQ], [1e-3], A.IrregularZenithAngleGrid(np.linspace(0.0, 180.0, 181))
)
two_series, _ = two_temperatures.to_tro_spectral_with_report([250.0, 300.0], [FREQ], 8)
two_species = A.ScatteringHabit(two_series, A.MonodispersePSD(drops), 1.0, 3.0)


def extinction(temperature):
    p = A.AtmPoint()
    p.temperature = temperature
    p.pressure = 1e5
    p[drops] = 1.0
    return np.asarray(two_species.get_bulk_scattering_properties_tro_spectral(p, [FREQ], 8).extinction_matrix)[0, 0]


k250, k300 = extinction(250.0), extinction(300.0)
for temperature in [220.0, 330.0]:
    linear = k250 + (temperature - 250.0) / 50.0 * (k300 - k250)
    assert abs(extinction(temperature) - linear) <= 1e-12 * abs(linear), (temperature, extinction(temperature), linear)
