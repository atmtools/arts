"""The spectral forms of scattering data through pyarts3: Legendre series and azimuthal Fourier modes.

- Gridded TRO data convert to a Legendre series with a report on how well it
  represents them; the series evaluates exactly on any grid.
- A habit of gridded data refuses to give a Legendre series, and names the
  conversion; converted, it gives one of the degree asked for.
- The azimuthal Fourier modes of the laboratory-frame phase matrix
  (get_bulk_scattering_properties_aro_fourier) evaluate to the gridded
  laboratory-frame phase matrix at any azimuth, exactly for Rayleigh
  scattering, whose modes above m = 2 vanish.
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
except RuntimeError as e:
    assert "to_tro_spectral" in str(e), str(e)
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
