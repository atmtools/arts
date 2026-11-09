"""Monte Carlo integration with sensor weighting and metadata."""

import numpy as np
import pyarts3 as pyarts

A = pyarts.arts
ws = pyarts.Workspace()

ws.surf_fieldPlanet(option="Earth")
ws.surf_field[A.SurfaceKey("t")] = 280.0
ws.atm_fieldRead(toa=30e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1)
ws.abs_speciesSet(species=[])
ws.abs_bands = {}
ws.scat_species = []
ws.spectral_propmat_agendaAuto()
ws.ray_path_observer_agendaSetGeometric(
    max_step_option="step", add_crossings=1, remove_non_atm=1
)
ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_rad_surface_agendaSet(option="Blackbody")


@pyarts.arts_agenda(ws=ws, fix=False)
def spectral_rad_observer_agenda(ws):
    ws.spectral_radMonteCarlo(mc_seed=9182, mc_min_iter=8, mc_max_iter=8)


freq = A.AscendingGrid([220e9, 230e9])
pos = [20e3, 0.0, 0.0]
ws.measurement_sensorInit()
ws.measurement_sensorAddSimple(freq_grid=freq, pos=pos, los=[180.0, 0.0])
ws.measurement_sensorAddSimple(freq_grid=freq, pos=pos, los=[0.0, 0.0])

surface_temperature_step = 0.5
ws.jac_targetsAddSurface(target="t", d=surface_temperature_step)
ws.jac_targetsFinalize()
ws.measurement_vecFromSensor()

expected_surface = [A.physics.planck(f, 280.0) for f in freq]
expected_space = [A.physics.planck(f, 2.725) for f in freq]
assert np.allclose(np.asarray(ws.measurement_vec), expected_surface + expected_space)
expected_jac = [
    (A.physics.planck(f, 280.0 + surface_temperature_step) - A.physics.planck(f, 280.0))
    / surface_temperature_step
    for f in freq
] + [0.0, 0.0]
assert np.allclose(np.asarray(ws.measurement_jac)[:, 0], expected_jac)
assert len(ws.measurement_sensor) == 2
assert all(len(sensor.elements) == 2 for sensor in ws.measurement_sensor)

surface, space = ws.measurement_sensor.split(ws.measurement_vec)
assert np.allclose(np.asarray(ws.measurement_sensor[0].reshape(surface).data), expected_surface)
assert np.allclose(np.asarray(ws.measurement_sensor[1].reshape(space).data), expected_space)

# A scattering-property target is perturbed through AtmField in model-state
# space.  Common random numbers make repeating the stochastic Jacobian exactly
# reproducible for a fixed seed.
ext = A.ScatteringSpeciesProperty("cloud", A.ParticulateProperty.Extinction)
ssa = A.ScatteringSpeciesProperty(
    "cloud", A.ParticulateProperty.SingleScatteringAlbedo
)
ws.atm_field[ext] = 1e-4
ws.atm_field[ssa] = 0.5
ws.scat_species = [A.HenyeyGreensteinScatterer(ext, ssa, 0.2)]
ws.jac_targetsInit()
ws.jac_targetsAddAtmosphere(target=ext, d=1e-5)
ws.jac_targetsAddSensorFrequencyPolyOffset(d=1e6, sensor_elem=0)
ws.jac_targetsFinalize()

ws.measurement_vecFromSensor()
ws.measurement_jacAddSensorJacobianPerturbations()
scattering_jac = np.asarray(ws.measurement_jac).copy()
assert np.all(np.isfinite(scattering_jac))
assert np.any(scattering_jac != 0.0)
assert np.any(scattering_jac[:2, 1] != 0.0)
assert np.all(scattering_jac[2:, 1] == 0.0)
ws.measurement_vecFromSensor()
ws.measurement_jacAddSensorJacobianPerturbations()
assert np.array_equal(scattering_jac, np.asarray(ws.measurement_jac))
