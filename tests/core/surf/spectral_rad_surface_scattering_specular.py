"""Plumbing tests for spectral_radSurfaceScatteringSpecular workspace method.

Specular-capable models (Fresnel, FlatScalar) exist, but these tests
deliberately register only a pure diffuse (Lambertian) model, so they
verify the plumbing:

1. Basic execution (smoke test): method runs and produces finite output
2. Zero-specular consistency: with a pure diffuse (Lambertian) model the
   specular BRDF and emissivity are zero, so the output radiance is zero
3. Jacobian shape: correct dimensions when jac_targets is non-empty
"""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts


def setup_workspace_base(freq_grid):
    """Create and configure a minimal workspace for the specular test."""
    ws = pyarts.Workspace()

    ws.freq_grid = freq_grid

    # Minimal atmosphere
    ws.abs_species = []
    ws.abs_bands = {}
    ws.spectral_propmat_agendaAuto()

    ws.atm_fieldInit(toa=100e3)

    # Surface field setup
    ws.surf_fieldEarth()
    ws.surf_field["t"] = 280.0

    # Jacobian targets (empty for basic tests)
    ws.jac_targets = arts.JacobianTargets()

    # Ray point at surface, nadir-looking (downward)
    ws.ray_point = arts.PropagationPathPoint()
    ws.ray_point.pos = [0.0, 0.0, 0.0]
    ws.ray_point.los = [180.0, 0.0]

    # suns is a required input of the scattering methods (the sun-beam exclusion
    # gate); empty list keeps the gate inert
    ws.suns = []

    # Agendas
    ws.spectral_rad_incoming_agendaSet(option="Emission")
    ws.ray_path_observer_agendaSetGeometric()

    return ws


def create_surface_models(freq_grid, reflectivity):
    """Create MapOfSurfaceScatteringModel with LambertianSurfaceScatterer."""
    refl_data = np.full(len(freq_grid), reflectivity)

    refl_field = arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[freq_grid],
        data=refl_data.tolist()
    )

    scatterer = arts.LambertianSurfaceScatterer(refl_field)

    surface_models = arts.MapOfSurfaceScatteringModel()
    surface_models.add("lambertian", scatterer)

    return surface_models


def add_surface_mask(ws, tag_key="lambertian"):
    """Add the surface property tag mask to ws.surf_field."""
    ws.surf_field[arts.SurfacePropertyTag(tag_key)] = 1.0


# ============================================================================
# Test 1: Smoke test / basic execution
# ============================================================================
def test_spectral_rad_surface_scattering_specular_basic():
    """Test basic execution with a diffuse (Lambertian) model."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringSpecular()

    assert len(ws.spectral_rad) == len(freq_grid), \
        f"spectral_rad size mismatch: {len(ws.spectral_rad)} != {len(freq_grid)}"

    rad_array = np.array([
        float(ws.spectral_rad[i][s]) for i in range(len(freq_grid)) for s in range(4)
    ])
    assert np.all(np.isfinite(rad_array)), \
        f"spectral_rad contains non-finite values: {rad_array}"

    print("Test 1 passed: basic execution with finite output")


# ============================================================================
# Test 2: Zero-specular consistency
# ============================================================================
def test_spectral_rad_surface_scattering_specular_zero_with_diffuse_model():
    """Lambertian models have zero specular BRDF and zero specular emissivity,
    so the specular method must produce zero radiance."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringSpecular()

    rad_array = np.array([float(ws.spectral_rad[i][0]) for i in range(len(freq_grid))])
    assert np.all(np.isclose(rad_array, 0.0)), \
        f"Specular radiance should be zero for a pure diffuse model: {rad_array}"

    print("Test 2 passed: pure diffuse model yields zero specular radiance")


# ============================================================================
# Test 3: Jacobian shape check
# ============================================================================
def test_spectral_rad_surface_scattering_specular_jacobian():
    """Test Jacobian computation with non-empty jac_targets."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    # Set up a Jacobian target: perturb the surface temperature.
    ws.measurement_sensorInit()
    ws.jac_targetsAddSurface(target="t")
    ws.jac_targetsFinalize()
    x_size = ws.jac_targets.x_size()

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringSpecular()

    jac_array = np.array(ws.spectral_rad_jac)
    assert jac_array.shape[0] == x_size, \
        f"Jacobian first dimension mismatch: {jac_array.shape[0]} != {x_size}"
    assert jac_array.shape[1] == len(freq_grid), \
        f"Jacobian second dimension mismatch: {jac_array.shape[1]} != {len(freq_grid)}"

    assert np.all(np.isfinite(jac_array)), \
        "Jacobian contains non-finite values"

    print("Test 3 passed: Jacobian has correct shape and finite values")


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    test_spectral_rad_surface_scattering_specular_basic()
    test_spectral_rad_surface_scattering_specular_zero_with_diffuse_model()
    test_spectral_rad_surface_scattering_specular_jacobian()
    print("\nAll tests passed!")
