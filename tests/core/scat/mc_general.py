"""Regression tests for the ARTS3-native passive Monte Carlo solver."""

import numpy as np
import pyarts3 as pyarts

A = pyarts.arts
ws = pyarts.Workspace()
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[A.SurfaceKey("t")] = 280.0
ws.atm_fieldRead(toa=30e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1)
ws.abs_speciesSet(species=[])
ws.spectral_propmat_agendaAuto()
ws.ray_path_observer_agendaSetGeometric(
    max_step_option="step", add_crossings=1, remove_non_atm=1
)
ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_rad_surface_agendaSet(option="Blackbody")
ws.mc_antenna.set_pencil_beam()

common = dict(
    frequency=230e9,
    sensor_pos=[20e3, 0.0, 0.0],
    sensor_los=[180.0, 0.0],
    mc_seed=1729,
    mc_min_iter=32,
    mc_max_iter=32,
    mc_max_scatorder=4,
)

# With no extinction every history reaches the same blackbody boundary.
ws.scat_species = []
ws.MCGeneral(**common)
surface = np.asarray(ws.mc_spectral_rad)
assert ws.mc_iteration_count == 32
assert np.allclose(np.asarray(ws.mc_error), np.zeros(4), atol=1e-22)
assert np.isclose(surface[0], A.physics.planck(common["frequency"], 280.0))
assert np.array_equal(surface[1:], np.zeros(3))

# Exercise collision, scattering, uncertainty, and deterministic seed behavior.
ext = A.ScatteringSpeciesProperty("cloud", A.ParticulateProperty.Extinction)
ssa = A.ScatteringSpeciesProperty(
    "cloud", A.ParticulateProperty.SingleScatteringAlbedo
)
grids = ws.atm_field["t"].data.grids
alt = np.asarray(grids[0])
cloud = (alt >= 5e3) & (alt <= 15e3)
ws.atm_field[ext] = A.GriddedField3(
    data=np.where(cloud, 2e-4, 0.0)[:, None, None], grids=grids
)
ws.atm_field[ssa] = A.GriddedField3(
    data=np.where(cloud, 0.8, 0.0)[:, None, None], grids=grids
)
ws.scat_species = [A.HenyeyGreensteinScatterer(ext, ssa, 0.2)]

scattering = common | {"mc_min_iter": 128, "mc_max_iter": 128}
ws.MCGeneral(**scattering)
rad = np.asarray(ws.mc_spectral_rad).copy()
err = np.asarray(ws.mc_error).copy()
assert np.all(np.isfinite(rad)) and np.all(np.isfinite(err))
assert rad[0] > 0.0 and err[0] > 0.0
# An unpolarizing HG scatterer must preserve an unpolarized boundary source;
# this still exercises the four-Stokes Mueller transport and basis transforms.
assert np.allclose(rad[1:], np.zeros(3), atol=1e-300)
ws.MCGeneral(**scattering)
assert np.array_equal(rad, np.asarray(ws.mc_spectral_rad))
assert np.array_equal(err, np.asarray(ws.mc_error))

# An absolute error threshold can stop after the requested minimum histories.
ws.MCGeneral(**(common | {"mc_min_iter": 10, "mc_max_iter": 100, "mc_std_err": 1.0}))
assert ws.mc_iteration_count == 10

# A homogeneous, optically thick atmosphere makes every history terminate in
# the same absorption source.  Exercise the additive, polarized NLTE emission
# supplied by spectral_propmat_agenda rather than reconstructing LTE in MCGeneral.
alpha = 1e-2
nlte_emission = np.array([2e-13, -3e-14, 4e-14, -5e-14])
temperature = 250.0
ws.atm_field["t"].data.data[:] = temperature
ws.scat_species = []


def constant_nlte(freq, wind_jac, jac, species, path_point, atm):
    local = pyarts.Workspace()
    local.spectral_propmatInit(freq_grid=freq, jac_targets=jac)
    local.spectral_propmat.K[:, 0] = alpha
    local.spectral_propmat.nlte[:] = nlte_emission
    local.spectral_propmat.nlte_on = True
    return (local.spectral_propmat,)


ws.spectral_propmat_agendaSetOperator(f=constant_nlte)
ws.MCGeneral(**(common | {"mc_min_iter": 16, "mc_max_iter": 16}))
expected = nlte_emission / alpha
expected[0] += A.physics.planck(common["frequency"], temperature)
assert np.allclose(np.asarray(ws.mc_spectral_rad), expected, rtol=1e-14)
assert np.all(np.asarray(ws.mc_error) <= 1e-8 * np.abs(expected))
