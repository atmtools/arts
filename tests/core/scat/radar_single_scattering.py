"""ARTS2-style deterministic active-radar forward-model regression tests."""

from pathlib import Path

import numpy as np
import pyarts3 as pyarts


A = pyarts.arts
ws = pyarts.Workspace()


def manual_atmospheric_derivative(workspace, key, flat_index, step, kwargs):
    """Central difference made only by explicit atmospheric-field edits."""
    original = np.asarray(workspace.atm_field[key].data.data).copy()
    workspace.jac_targetsOff()
    try:
        plus = original.copy()
        plus.flat[flat_index] += step
        workspace.atm_field[key].data.data = plus
        workspace.measurement_vecFromRadarSingleScattering(**kwargs)
        y_plus = np.asarray(workspace.measurement_vec).copy()

        minus = original.copy()
        minus.flat[flat_index] -= step
        workspace.atm_field[key].data.data = minus
        workspace.measurement_vecFromRadarSingleScattering(**kwargs)
        y_minus = np.asarray(workspace.measurement_vec).copy()
    finally:
        workspace.atm_field[key].data.data = original
    return (y_plus - y_minus) / (2.0 * step)

ws.surf_fieldPlanet(option="Earth")
ws.surf_field[A.SurfaceKey("h")] = 0.0
ws.surf_field[A.SurfaceKey("t")] = 280.0
ws.atm_fieldRead(
    toa=20e3,
    basename="planets/Earth/afgl/tropical/",
    missing_is_zero=1,
)
ws.abs_speciesSet(species=[])
ws.abs_bands = {}
ws.spectral_propmat_agendaAuto()
ws.max_stepsize = 250.0
ws.ray_path_observer_agendaSetGeometric(
    max_step_option="step", add_crossings=1, remove_non_atm=1
)

# A scalar, isotropic scatterer makes the native backscatter coefficient
# exactly extinction * single-scattering-albedo / (4 pi).  This checks the
# same normalization on which ARTS2 TestIyActive's -30 dBZe case depends.
extinction = A.ScatteringSpeciesProperty(
    "cloud", A.ParticulateProperty.Extinction
)
ssa = A.ScatteringSpeciesProperty(
    "cloud", A.ParticulateProperty.SingleScatteringAlbedo
)
grids = ws.atm_field["t"].data.grids
altitude = np.asarray(grids[0])
cloud = (altitude >= 5e3) & (altitude <= 15e3)
extinction_value = 1e-4
ssa_value = 0.8
ws.atm_field[extinction] = A.GriddedField3(
    data=np.where(cloud, extinction_value, 0.0)[:, None, None], grids=grids
)
ws.atm_field[ssa] = A.GriddedField3(
    data=np.where(cloud, ssa_value, 0.0)[:, None, None], grids=grids
)
ws.scat_species = [A.HenyeyGreensteinScatterer(extinction, ssa, 0.0)]

# ARTS2 polarization index 5 is [I,Q]/2 on receive and [I,Q] on transmit.
# The helper creates one standard SensorObsel per range bin, so the radar
# result can be used directly as an OEM measurement vector.
range_bins = np.arange(0.0, 20.001e3, 500.0)
ws.freq_grid = [94e9]
ws.measurement_sensorInit()
ws.radar_range_limits = np.empty((0, 2))
ws.measurement_sensorAddSimpleRadar(
    pos=[20e3, 0.0, 0.0],
    los=[180.0, 0.0],
    pol=[0.5, 0.5, 0.0, 0.0],
    range_bins=range_bins,
)

common = dict(
    transmitted_stokes=[1.0, 1.0, 0.0, 0.0],
    range_mode="Altitude",
    pext_scaling=0.0,
)
aux_names = [
    "Radiative background",
    "Backscattering",
    "Abs species extinction",
    "Particle extinction",
]
ws.measurement_vecFromRadarSingleScattering(
    **common, unit="1", aux_vars=aux_names
)
native = np.asarray(ws.measurement_vec).copy()
aux = np.asarray(ws.radar_aux).copy()

expected_native = extinction_value * ssa_value / (4.0 * np.pi)
interior = (range_bins[:-1] >= 5.5e3) & (range_bins[1:] <= 14.5e3)
assert np.allclose(native[interior], expected_native, rtol=2e-13)
assert np.allclose(aux[1, interior], expected_native, rtol=2e-13)
assert np.allclose(aux[2, interior], 0.0)
# The extinction auxiliary is independent of the 0.5 I coefficient in the
# polarization selector, matching old yRadar behavior.
assert np.allclose(aux[3, interior], 0.0)  # pext_scaling is zero

# The round-trip time of every gate of this vertical vacuum path.
centres = 0.5 * (range_bins[1:] + range_bins[:-1])
np.testing.assert_allclose(
    ws.radar_time, 2.0 * (20e3 - centres) / A.constants.c, rtol=1e-12
)

# The Magnus operator differs only for non-commuting propagation matrices.
ws.rte_option = "magop"
ws.measurement_vecFromRadarSingleScattering(**common, unit="1")
np.testing.assert_allclose(ws.measurement_vec, native, rtol=1e-12, atol=1e-30)
ws.rte_option = "linsrc"

# Ze and dBZe conversions, including the old clipping convention.
k2 = 0.93
ws.measurement_vecFromRadarSingleScattering(**common, unit="Ze", k2=k2)
ze = np.asarray(ws.measurement_vec).copy()
wavelength = A.constants.c / ws.freq_grid[0]
ze_factor = 4e18 * wavelength**4 / (np.pi**4 * k2)
assert np.allclose(ze[interior], ze_factor * native[interior], rtol=2e-13)

ws.measurement_vecFromRadarSingleScattering(
    **common, unit="dBZe", k2=k2, dbze_min=-80.0
)
dbze = np.asarray(ws.measurement_vec).copy()
assert np.allclose(dbze[interior], 10.0 * np.log10(ze[interior]))
assert np.all(dbze[native == 0.0] == -80.0)

# Distance and round-trip-time gates must select the same cloud as altitude
# gates for this vertical vacuum path.  Reverse the altitude edges because
# one-way distance increases downwards from the radar.
distance_bins = 20e3 - range_bins[::-1]
ws.radar_range_limits = np.column_stack(
    (distance_bins[:-1], distance_bins[1:])
)
ws.measurement_vecFromRadarSingleScattering(
    **(common | {"range_mode": "Distance"}), unit="1"
)
distance_signal = np.asarray(ws.measurement_vec).copy()
assert np.allclose(distance_signal, native[::-1], rtol=2e-12, atol=1e-30)

time_bins = 2.0 * distance_bins / A.constants.c
ws.radar_range_limits = np.column_stack((time_bins[:-1], time_bins[1:]))
ws.measurement_vecFromRadarSingleScattering(
    **(common | {"range_mode": "RoundTripTime"}), unit="1"
)
assert np.allclose(
    np.asarray(ws.measurement_vec), distance_signal, rtol=2e-12, atol=1e-30
)

# Restore altitude gates and compare a mapped scattering-property Jacobian
# against an explicit central difference.  This is the ARTS3 equivalent of
# TestIyActive_wfuns' particle-profile Jacobian check.
ws.radar_range_limits = np.column_stack((range_bins[:-1], range_bins[1:]))
step = 1e-6
ws.jac_targetsInit()
ws.jac_targetsAddAtmosphere(target=extinction, d=step)
ws.jac_targetsFinalize()
ws.measurement_vecFromRadarSingleScattering(**common, unit="1")
jac = np.asarray(ws.measurement_jac).copy()
assert jac.shape == (len(ws.measurement_sensor), altitude.size)
assert np.all(np.isfinite(jac))

iz = int(np.argmin(np.abs(altitude - 10e3)))
base_extinction = np.asarray(ws.atm_field[extinction].data.data).copy()
ws.jac_targetsOff()

plus = base_extinction.copy()
plus[iz, 0, 0] += step
ws.atm_field[extinction].data.data = plus
ws.measurement_vecFromRadarSingleScattering(**common, unit="1")
y_plus = np.asarray(ws.measurement_vec).copy()

minus = base_extinction.copy()
minus[iz, 0, 0] -= step
ws.atm_field[extinction].data.data = minus
ws.measurement_vecFromRadarSingleScattering(**common, unit="1")
y_minus = np.asarray(ws.measurement_vec).copy()

explicit = (y_plus - y_minus) / (2.0 * step)
assert np.allclose(jac[:, iz], explicit, rtol=2e-10, atol=1e-15)


# Run the central TestIyActive normalization check with an in-memory Mie
# droplet matching the original ARTS2 material model and particle parameters.
data = Path(__file__).parent / "data" / "radar_arts2"
old = pyarts.Workspace()
old.abs_speciesSet(
    species=["N2-SelfContStandardType", "O2-PWR98", "H2O-PWR98"]
)
old.abs_bands = {}
old.surf_fieldPlanet(option="Earth")
old.surf_field[A.SurfaceKey("h")] = 0.0
old.surf_field[A.SurfaceKey("t")] = 273.15
old.atm_fieldRead(
    toa=16e3,
    basename="planets/Earth/afgl/tropical/",
    missing_is_zero=1,
)
old.atm_field["t"].data.data[:] = 273.15

number_density = A.ScatteringSpeciesProperty(
    "arts2_droplet", A.ParticulateProperty.NumberDensity
)
fine_altitude = np.arange(0.0, 16.001e3, 50.0)
old.atm_field[number_density] = A.GriddedField3(
    data=np.where(
        (fine_altitude >= 50.0) & (fine_altitude <= 5e3), 6.4e4, 0.0
    )[:, None, None],
    grids=[fine_altitude, *old.atm_field["t"].data.grids[1:]],
)

# Same Liebe (1993) water model as Atmlab eps_water_liebe93.m. Its
# relaxation-frequency coefficient is 146, not epswater93.m's 146.4.
# Preserve the original frequency, temperatures, diameter and angular grid.
temperatures = np.array([268.15, 273.15, 278.15])
theta = 300.0 / temperatures - 1.0
e0 = 77.66 + 103.3 * theta
e1 = 0.0671 * e0
g1 = 20.2 - 146.0 * theta + 316.0 * theta**2
permittivity = e0 - 94.0 * (
    (e0 - e1) / (94.0 + 1j * g1)
    + (e1 - 3.52) / (94.0 + 1j * 39.8 * g1)
)
habit = A.ParticleHabit.sphere(
    temperatures, [94e9], [50e-6],
    A.IrregularZenithAngleGrid(np.arange(181.0)),
    np.sqrt(permittivity)[:, None], 1000.0,
)
particle = habit[0]
np.testing.assert_allclose(particle.properties.mass, 1000 * np.pi / 6 * (50e-6)**3)
phase = np.asarray(particle.phase_matrix)
np.testing.assert_array_equal(particle.backscatter_matrix, phase[:, :, -1, :])
np.testing.assert_array_equal(particle.forwardscatter_matrix, phase[:, :, 0, :])
# Endpoint data must also be exact when the sampled phase grid omits 0/180.
sparse = A.ParticleHabit.sphere(
    temperatures, [94e9], [50e-6],
    A.IrregularZenithAngleGrid([30.0, 90.0, 150.0]),
    np.sqrt(permittivity)[:, None], 1000.0,
)[0]
np.testing.assert_array_equal(sparse.backscatter_matrix, particle.backscatter_matrix)
np.testing.assert_array_equal(sparse.forwardscatter_matrix, particle.forwardscatter_matrix)
liquid = A.ParticleHabit.liquid_sphere(
    temperatures, [94e9], [50e-6], A.IrregularZenithAngleGrid(np.arange(181.0))
)[0]
np.testing.assert_allclose(liquid.properties.mass, particle.properties.mass)
np.testing.assert_array_equal(liquid.backscatter_matrix, np.asarray(liquid.phase_matrix)[:, :, -1, :])
np.testing.assert_array_equal(liquid.forwardscatter_matrix, np.asarray(liquid.phase_matrix)[:, :, 0, :])
old.scat_species = [
    A.ScatteringHabit(
        habit,
        A.MonodispersePSD(number_density, 0.0, 400.0),
        1000.0 * np.pi / 6,
        3.0,
    )
]

old.max_stepsize = 50.0
old.ray_path_observer_agendaSetGeometric(
    max_step_option="step", add_crossings=1, remove_non_atm=1
)
old.freq_grid = [94e9]
old.measurement_sensorInit()
old.radar_range_limits = np.empty((0, 2))
old.measurement_sensorAddSimpleRadar(
    pos=[100e3, 0.0, 0.0],
    los=[180.0, 0.0],
    pol=[0.5, 0.5, 0.0, 0.0],
    range_bins=np.arange(0.0, 10.001e3, 500.0),
)

old_common = dict(
    transmitted_stokes=[1.0, 1.0, 0.0, 0.0],
    range_mode="Altitude",
    ze_tref=273.15,
    dbze_min=-80.0,
)

# With extinction disabled the old test expects a -30 dBZe peak to 0.005 dB.
old.abs_speciesSet(species=[])
old.spectral_propmat_agendaAuto()
old.measurement_vecFromRadarSingleScattering(
    **old_common, unit="dBZe", pext_scaling=0.0
)
assert abs(np.nanmax(old.measurement_vec) + 30.0) < 0.005

# Particle extinction retains the old 0.01 dB peak tolerance.
old.measurement_vecFromRadarSingleScattering(
    **old_common, unit="dBZe", pext_scaling=1.0
)
particle_peak = np.nanmax(old.measurement_vec)
assert abs(particle_peak + 30.0) < 0.01

# The scattering-habit/PSD derivative is analytical in the forward method.
# Build its numerical reference here by manually perturbing one number-density
# grid point and rerunning with Jacobians disabled.
particle_step = 10.0
particle_index = int(np.argmin(np.abs(fine_altitude - 2.5e3)))
old.jac_targetsInit()
old.jac_targetsAddAtmosphere(target=number_density, d=particle_step)
old.jac_targetsFinalize()
old.measurement_vecFromRadarSingleScattering(
    **old_common, unit="Ze", pext_scaling=1.0
)
particle_jac = np.asarray(old.measurement_jac).copy()
particle_explicit = manual_atmospheric_derivative(
    old,
    number_density,
    particle_index,
    particle_step,
    old_common | {"unit": "Ze", "pext_scaling": 1.0},
)
assert np.allclose(
    particle_jac[:, particle_index],
    particle_explicit,
    rtol=2e-7,
    atol=1e-10,
)

# Restore the old predefined gases and verify that gas attenuation and all
# four legacy auxiliary outputs participate in the calculation.
old.abs_speciesSet(
    species=["N2-SelfContStandardType", "O2-PWR98", "H2O-PWR98"]
)
old.abs_predef_dataReadSpeciesSplitCatalog(
    basename=str(data) + "/", name_missing=1, ignore_missing=1
)
old.spectral_propmat_agendaAuto()
old.measurement_vecFromRadarSingleScattering(
    **old_common,
    unit="dBZe",
    pext_scaling=1.0,
    aux_vars=aux_names,
)
assert np.nanmax(old.measurement_vec) < particle_peak
old_aux = np.asarray(old.radar_aux)
assert old_aux.shape == (4, len(old.measurement_sensor))
assert np.nanmax(old_aux[2]) > 0.0
assert np.nanmax(old_aux[3]) > 0.0

# TestIyActive_wfuns also exercises absorbing-species and temperature
# derivatives.  ARTS3 assembles these through the finalized atmospheric state
# mapping, so both blocks must be finite and nonzero in standard OEM layout.
old.jac_targetsInit()
old.jac_targetsAddSpeciesVMR(species="H2O", d=1e-5)
old.jac_targetsAddTemperature(d=0.1)
old.jac_targetsFinalize()
old.measurement_vecFromRadarSingleScattering(
    **old_common, unit="Ze", pext_scaling=1.0
)
old_jac = np.asarray(old.measurement_jac).copy()
nh2o = np.asarray(old.atm_field["H2O"].data.data).size
assert old_jac.shape[0] == len(old.measurement_sensor)
assert np.all(np.isfinite(old_jac))
assert np.max(np.abs(old_jac[:, :nh2o])) > 0.0
assert np.max(np.abs(old_jac[:, nh2o:])) > 0.0

# Independently validate both propagation-matrix derivative blocks.  These
# finite differences deliberately live in the test: the production radar
# method receives analytical derivatives from spectral_propmat_agenda and
# propagates them through the two-way rtepack transmission.
h2o_altitude = np.asarray(old.atm_field["H2O"].data.grids[0])
h2o_index = int(np.argmin(np.abs(h2o_altitude - 2.5e3)))
gas_kwargs = old_common | {"unit": "Ze", "pext_scaling": 1.0}
h2o_explicit = manual_atmospheric_derivative(
    old, "H2O", h2o_index, 1e-7, gas_kwargs
)
temperature_explicit = manual_atmospheric_derivative(
    old, "t", h2o_index, 0.05, gas_kwargs
)
assert np.allclose(
    old_jac[:, h2o_index], h2o_explicit, rtol=3e-4, atol=1e-10
)
assert np.allclose(
    old_jac[:, nh2o + h2o_index],
    temperature_explicit,
    rtol=2e-3,
    atol=1e-10,
)
