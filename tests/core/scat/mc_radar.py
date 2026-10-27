"""Seeded regression tests for the ARTS3-native MCRadar interface."""

import numpy as np
import pyarts3 as pyarts

A = pyarts.arts

ws = pyarts.Workspace()
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[A.SurfaceKey("t")] = 280.0
ws.atm_fieldRead(
    toa=30e3,
    basename="planets/Earth/afgl/tropical/",
    missing_is_zero=1,
)
ws.abs_speciesSet(species=[])
ws.spectral_propmat_agendaAuto()
ws.ray_path_observer_agendaSetGeometric(
    max_step_option="step",
    add_crossings=1,
    remove_non_atm=1,
)

# A horizontally homogeneous, isotropic cloud.  Extinction is deliberately
# small enough that both the cloud top and bottom remain visible.
extinction = A.ScatteringSpeciesProperty(
    "cloud", A.ParticulateProperty.Extinction
)
ssa = A.ScatteringSpeciesProperty(
    "cloud", A.ParticulateProperty.SingleScatteringAlbedo
)
grids = ws.atm_field["t"].data.grids
altitude = np.asarray(grids[0])
cloud = (altitude >= 5e3) & (altitude <= 15e3)
ws.atm_field[extinction] = A.GriddedField3(
    data=np.where(cloud, 1e-4, 0.0)[:, None, None], grids=grids
)
ws.atm_field[ssa] = A.GriddedField3(
    data=np.where(cloud, 0.8, 0.0)[:, None, None], grids=grids
)
ws.scat_species = [A.HenyeyGreensteinScatterer(extinction, ssa, 0.0)]

ws.mc_antenna.set_gaussian_fwhm(0.6, 0.6)
kwargs = dict(
    frequency=94e9,
    sensor_pos=[20e3, 0.0, 0.0],
    sensor_los=[180.0, 0.0],
    mc_y_tx=[1.0, 1.0, 0.0, 0.0],
    range_bins=np.arange(0.0, 20.001e3, 1e3),
    mc_seed=1729,
    mc_max_iter=64,
    mc_max_scatorder=1,
    unit="1",
)

ws.MCRadar(**kwargs)
signal = np.array(ws.radar_signal)
error = np.array(ws.radar_error)
gaussian_time = np.array(ws.radar_time)

assert isinstance(ws.radar_signal, A.StokvecVector)
assert isinstance(ws.radar_error, A.StokvecVector)
assert signal.shape == (20, 4)
assert error.shape == signal.shape
assert np.all(np.isfinite(signal))
assert np.all(np.isfinite(error))
assert np.all(error >= 0.0)
assert np.any(signal[:, 0] > 0.0)
assert gaussian_time.shape == (20,)
assert np.all(np.diff(gaussian_time) > 0.0)

# The same seed must give bitwise-identical Monte Carlo output.
ws.MCRadar(**kwargs)
assert np.array_equal(signal, np.asarray(ws.radar_signal))
assert np.array_equal(error, np.asarray(ws.radar_error))

# Ze is a linear conversion of the native return, including its uncertainty.
ws.MCRadar(**(kwargs | {"unit": "Ze"}))
ze = np.array(ws.radar_signal)
ze_error = np.array(ws.radar_error)
mask = signal != 0.0
factor = ze[mask] / signal[mask]
assert np.allclose(factor, factor[0], rtol=1e-13)
error_mask = error != 0.0
assert np.allclose(ze_error[error_mask] / error[error_mask], factor[0])

# In the pencil-beam, single-scattering limit, MCRadar must reproduce the
# ARTS3 deterministic rtepack radar result.  One history is sufficient because
# there is then no random angular sampling.  Both methods share the gate
# quadrature and the two-way transmission, so they agree to rounding, also in
# Ze.
ws.mc_antenna.set_pencil_beam()
ws.max_stepsize = 100.0
ws.MCRadar(**(kwargs | {"mc_max_iter": 1}))
mc_pencil = np.asarray(ws.radar_signal)[:, 0].copy()
mc_time = np.asarray(ws.radar_time).copy()
ws.MCRadar(**(kwargs | {"mc_max_iter": 1, "unit": "Ze"}))
mc_pencil_ze = np.asarray(ws.radar_signal)[:, 0].copy()

ws.freq_grid = [kwargs["frequency"]]
ws.measurement_sensorInit()
ws.radar_range_limits = np.empty((0, 2))
ws.measurement_sensorAddSimpleRadar(
    pos=kwargs["sensor_pos"],
    los=kwargs["sensor_los"],
    pol=[1.0, 0.0, 0.0, 0.0],
    range_bins=kwargs["range_bins"],
)
ws.jac_targetsOff()
deterministic_kwargs = dict(
    transmitted_stokes=kwargs["mc_y_tx"],
    range_mode="Distance",
    pext_scaling=1.0,
    k2=0.93,
)
ws.measurement_vecFromRadarSingleScattering(**deterministic_kwargs, unit="1")
deterministic = np.asarray(ws.measurement_vec).copy()
np.testing.assert_allclose(mc_pencil, deterministic, rtol=1e-10, atol=1e-30)
np.testing.assert_allclose(mc_time, ws.radar_time, rtol=1e-12)

ws.measurement_vecFromRadarSingleScattering(**deterministic_kwargs, unit="Ze")
np.testing.assert_allclose(mc_pencil_ze, ws.measurement_vec, rtol=1e-10, atol=1e-30)

# The round-trip time of a vertical path in vacuum is twice the distance to
# the gate centre divided by the speed of light.
centres = 0.5 * (kwargs["range_bins"][1:] + kwargs["range_bins"][:-1])
np.testing.assert_allclose(mc_time, 2.0 * centres / A.constants.c, rtol=1e-12)

# A narrow beam is uniformly filled by this horizontally homogeneous cloud.
# The antenna-weighted mean must then equal the pencil-beam return.
ws.mc_antenna.set_gaussian_fwhm(0.6, 0.6)
ws.MCRadar(**(kwargs | {"mc_max_iter": 16}))
gaussian = np.asarray(ws.radar_signal)[:, 0]
inside = (centres > 6e3) & (centres < 14e3)
np.testing.assert_allclose(gaussian[inside], mc_pencil[inside], rtol=1e-4)

# Invalid higher scattering orders are rejected explicitly until the collision
# sampler is added, rather than silently returning a single-scattering result.
try:
    ws.MCRadar(**(kwargs | {"mc_max_scatorder": 2}))
except RuntimeError:
    pass
else:
    raise AssertionError("MCRadar accepted an unsupported scattering order")
