"""Generated ARO radar onion-peeling and retrieval-state regression test."""

import numpy as np
import pyarts3 as pyarts


A = pyarts.arts
ws = pyarts.Workspace()

# Build the complete atmosphere in memory.  In particular, this test has no
# scattering XML, inversion table, catalog, or dependency on tmp/.
altitude = np.arange(0.0, 12.001e3, 1e3)
grids = [altitude, [0.0], [0.0]]
names = ["altitude", "latitude", "longitude"]


def field(profile):
    return A.GeodeticField3(
        data=np.asarray(profile)[:, None, None], grids=grids, grid_names=names
    )


ws.atm_field.top_of_atmosphere = altitude[-1]
ws.atm_field["t"] = field(260.0 + 20.0 * (1.0 - altitude / altitude[-1]))
ws.atm_field["p"] = field(1e5 * np.exp(-altitude / 7e3))
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[A.SurfaceKey("h")] = 0.0
ws.surf_field[A.SurfaceKey("t")] = 280.0
ws.abs_speciesSet(species=[])
ws.abs_bands = {}
ws.spectral_propmat_agendaAuto()
ws.max_stepsize = 100.0
ws.ray_path_observer_agendaSetGeometric(
    max_step_option="step", add_crossings=1, remove_non_atm=1
)

extinction = A.ScatteringSpeciesProperty(
    "generated-cloud", A.ParticulateProperty.Extinction
)
ssa = A.ScatteringSpeciesProperty(
    "generated-cloud", A.ParticulateProperty.SingleScatteringAlbedo
)

# HG is evaluated through its laboratory-frame (ARO) bulk properties by the
# radar model.  The varying profile gives both local backscatter and two-way
# particulate attenuation something nontrivial to retrieve.
truth = 2.5e-5 + 1.2e-5 * np.sin(np.pi * altitude / altitude[-1]) ** 2
truth[[0, -1]] = 0.0
ws.atm_field[extinction] = field(truth)
ws.atm_field[ssa] = field(np.full(altitude.size, 0.82))
ws.scat_species = [A.HenyeyGreensteinScatterer(extinction, ssa, 0.15)]

ws.freq_grid = [35e9]
ws.measurement_sensorInit()
ws.radar_range_limits = np.empty((0, 2))
range_edges = np.arange(0.5e3, 11.501e3, 1e3)
ws.measurement_sensorAddSimpleRadar(
    pos=[12e3, 0.0, 0.0],
    los=[180.0, 0.0],
    pol=[0.5, 0.5, 0.0, 0.0],
    range_bins=range_edges,
)

forward = dict(
    transmitted_stokes=[1.0, 1.0, 0.0, 0.0],
    range_mode="Altitude",
    unit="1",
    pext_scaling=1.0,
)
ws.jac_targetsOff()
ws.measurement_vecFromRadarSingleScattering(**forward)
observed = np.asarray(ws.measurement_vec).copy()
assert np.all(np.isfinite(observed)) and np.all(observed > 0.0)

# Supply a deliberately flat, low initial profile.  The boundary values are
# outside the range gates and remain an ordinary part of the state vector.
initial = np.full(altitude.size, 1.0e-5)
initial[[0, -1]] = 0.0
ws.atm_field[extinction] = field(initial)
ws.jac_targetsInit()
ws.jac_targetsAddAtmosphere(target=extinction, d=1e-7)
ws.jac_targetsFinalize()
ws.measurement_vec = observed

ws.model_state_vecFromRadarOnionPeeling(
    **forward,
    state_min=0.0,
    state_max=1e-3,
    max_step=2e-4,
    tolerance=2e-8,
    max_iterations=12,
    max_sweeps=12,
)

fit = np.asarray(ws.measurement_vec_fit)
state = np.asarray(ws.model_state_vec)
jac = np.asarray(ws.measurement_jac)
retrieved = np.asarray(ws.atm_field[extinction].data.data).ravel()

assert state.shape == truth.shape
assert jac.shape == (observed.size, truth.size)
assert np.all(np.isfinite(jac))
assert np.allclose(retrieved, state)
assert np.all(retrieved >= 0.0)
assert np.allclose(fit, observed, rtol=3e-4, atol=1e-13)

# The returned state is not an onion-specific object: applying it through the
# standard model-state method reproduces the returned field and radar fit.
saved = retrieved.copy()
ws.atm_field[extinction] = field(initial)
ws.atm_fieldFromModelState(model_state_targets=ws.jac_targets)
assert np.allclose(np.asarray(ws.atm_field[extinction].data.data).ravel(), saved)
ws.measurement_vecFromRadarSingleScattering(**forward)
assert np.allclose(np.asarray(ws.measurement_vec), fit, rtol=2e-12, atol=1e-15)
