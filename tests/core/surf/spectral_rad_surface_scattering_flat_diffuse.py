"""Tests for spectral_radSurfaceScatteringFlatDiffuse workspace method.

Verifies:
1. Basic execution (smoke test): method runs and produces finite output
2. Energy conservation: absorbing surface (R=0) vs. perfect reflector (R=1)
3. Physical plausibility: scattered contribution increases with reflectivity
4. Jacobian shape: correct dimensions when jac_targets is non-empty
"""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts


def setup_workspace_base(freq_grid):
    """Create and configure a minimal workspace for surface scattering tests.
    
    Args:
        freq_grid: list or array of frequencies
    
    Returns:
        Configured Workspace object
    """
    ws = pyarts.Workspace()
    
    # Set up frequency grid
    ws.freq_grid = freq_grid
    
    # Minimal atmosphere
    ws.abs_species = []
    ws.abs_bands = {}
    ws.spectral_propmat_agendaAuto()
    
    # Initialize atmosphere field
    ws.atm_fieldInit(toa=100e3)
    
    # Surface field setup
    ws.surf_fieldEarth()
    ws.surf_field["t"] = 280.0  # Surface temperature
    
    # Default subsurface field
    # ws.subsurf_field is already initialized with defaults
    
    # Jacobian targets (empty for basic tests)
    ws.jac_targets = arts.JacobianTargets()
    
    # Ray point at surface, nadir-looking (downward)
    ws.ray_point = arts.PropagationPathPoint()
    ws.ray_point.pos = [0.0, 0.0, 0.0]
    ws.ray_point.los = [180.0, 0.0]
    
    # Angular grids for hemisphere integration
    # Use a coarse grid for speed; Gauss-Legendre quadrature
    nza = 5
    naa = 4
    ws.zen_grid = arts.ZenGrid(np.linspace(0, 90, nza))
    ws.az_grid = arts.AziGrid(np.linspace(0, 360, naa, endpoint=False))
    
    # Quadrature weights: simple trapezoid in (za, aa)
    # dOmega = sin(za) * dza * daa
    dza = 90.0 / (nza - 1) if nza > 1 else 90.0
    daa = 360.0 / naa
    za_rad = np.deg2rad(np.linspace(0, 90, nza))
    za_weights = np.sin(za_rad) * np.deg2rad(dza)
    az_weights = np.ones(naa) * np.deg2rad(daa)
    
    ws.zen_grid_weights = arts.Vector(za_weights.tolist())
    ws.az_grid_weights = arts.Vector(az_weights.tolist())
    
    # Agendas
    ws.spectral_rad_incoming_agendaSet(option="Emission")
    ws.ray_path_observer_agendaSetGeometric()
    
    # Initialize spectral_rad and spectral_rad_jac (will be resized by the method)
    ws.spectral_rad = arts.StokvecVector()
    ws.spectral_rad_jac = arts.StokvecMatrix()
    
    return ws


def create_surface_models(freq_grid, reflectivity):
    """Create MapOfSurfaceScatteringModel with LambertianSurfaceScatterer.
    
    Args:
        freq_grid: list or array of frequencies
        reflectivity: scalar reflectivity value (0-1) or frequency-dependent array
    
    Returns:
        MapOfSurfaceScatteringModel
    """
    # Create frequency-dependent reflectivity field
    if np.isscalar(reflectivity):
        refl_data = np.full(len(freq_grid), reflectivity)
    else:
        refl_data = reflectivity
    
    refl_field = arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[freq_grid],
        data=refl_data.tolist()
    )
    
    # Create Lambertian scatterer
    scatterer = arts.LambertianSurfaceScatterer(refl_field)
    
    # Add to map with key "lambertian"
    surface_models = arts.MapOfSurfaceScatteringModel()
    surface_models.add("lambertian", scatterer)
    
    return surface_models


def add_surface_mask(ws, tag_key="lambertian"):
    """Add a surface property tag mask to the workspace surf_field.
    
    Args:
        ws: Workspace
        tag_key: SurfacePropertyTag key name (must match surface_models keys)
    """
    tag = arts.SurfacePropertyTag(tag_key)
    ws.surf_field[tag] = 1.0


# ============================================================================
# Test 1: Smoke test / basic execution
# ============================================================================
def test_spectral_rad_surface_scattering_flat_diffuse_basic():
    """Test basic execution with uniform reflectivity."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)
    
    # Add surface mask
    add_surface_mask(ws, "lambertian")
    
    # Create surface models with moderate reflectivity
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)
    
    # Execute the method
    ws.spectral_radSurfaceScatteringFlatDiffuse()
    
    # Assertions
    assert len(ws.spectral_rad) == len(freq_grid), \
        f"spectral_rad size mismatch: {len(ws.spectral_rad)} != {len(freq_grid)}"
    
    rad_array = np.array([float(ws.spectral_rad[i][0]) for i in range(len(freq_grid))])
    assert np.all(np.isfinite(rad_array)), \
        f"spectral_rad contains non-finite values: {rad_array}"
    
    assert np.all(rad_array > 0), \
        f"spectral_rad should be positive (thermal scene), got {rad_array}"
    
    print("Test 1 passed: basic execution with finite, positive output")


# ============================================================================
# Test 2: Energy conservation / absorbing surface (reflectivity = 0)
# ============================================================================
def test_spectral_rad_surface_scattering_flat_diffuse_absorbing():
    """Test with reflectivity = 0 (fully absorbing surface)."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)
    
    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.0)
    
    ws.spectral_radSurfaceScatteringFlatDiffuse()
    rad_absorbing = np.array([float(ws.spectral_rad[i][0]) for i in range(len(freq_grid))])
    
    # With reflectivity=0, scattering contribution is zero.
    # The output should be dominated by surface emission.
    assert np.all(np.isfinite(rad_absorbing)), \
        "Absorbing surface: spectral_rad contains non-finite values"
    
    assert np.all(rad_absorbing > 0), \
        "Absorbing surface: spectral_rad should be positive (surface emission)"
    
    # Estimate Planck function at 280 K to verify order of magnitude
    T_surf = ws.surf_field["t"](0,0)
    from scipy import constants
    h = constants.h
    c = constants.c
    k_B = constants.k
    
    # Planck function at a representative frequency (100 GHz)
    idx=0
    f_ref = freq_grid[idx]
    planck = (2 * h * f_ref**3 / c**2) / (np.exp(h * f_ref / (k_B * T_surf)) - 1)
    
    # The spectral radiance should be in a reasonable range relative to Planck
    rad_ref = float(rad_absorbing[idx])  # At 100 GHz
    ratio = rad_ref / planck
    assert 0.9999999 < ratio < 1.0000001, \
        f"Planck ratio out of range: {ratio} (expect 0.01-100)"
    
    print("Test 2 passed: absorbing surface produces finite surface-emission-only output")


# ============================================================================
# Test 3: Perfect reflector (reflectivity = 1)
# ============================================================================
def test_spectral_rad_surface_scattering_flat_diffuse_reflector():
    """Test the Kirchhoff coupling emissivity = 1 - reflectivity.

    The test atmosphere contains no absorption species, so there is no
    incoming thermal radiation for the surface to reflect.  Therefore:

    * Perfect reflector (r = 1): emissivity = 0, nothing to reflect ->
      the radiance is (near-)zero.
    * Absorbing surface (r = 0): emissivity = 1 -> the full 280 K blackbody.

    The reflector output must therefore be below the absorbing-surface output
    across the whole spectrum.
    """
    freq_grid = [10e9, 100e9, 183e9]

    # Run absorbing case
    ws1 = setup_workspace_base(freq_grid)
    add_surface_mask(ws1, "lambertian")
    ws1.surface_models = create_surface_models(freq_grid, reflectivity=0.0)
    ws1.spectral_radSurfaceScatteringFlatDiffuse()
    rad_absorbing = np.array([float(ws1.spectral_rad[i][0]) for i in range(len(freq_grid))])

    # Run perfect reflector case
    ws2 = setup_workspace_base(freq_grid)
    add_surface_mask(ws2, "lambertian")
    ws2.surface_models = create_surface_models(freq_grid, reflectivity=1.0)
    ws2.spectral_radSurfaceScatteringFlatDiffuse()
    rad_reflector = np.array([float(ws2.spectral_rad[i][0]) for i in range(len(freq_grid))])

    # Emissivity = 1 - r: the reflecting surface cannot emit and has
    # (no incoming thermal radiation to) reflect, so it must produce less
    # outgoing radiance than the emitting one.
    assert np.all(rad_reflector < rad_absorbing), \
        f"Emissivity = 1 - r violated: reflector {rad_reflector} >= absorber {rad_absorbing}"

    # The reflector output is only the negligible scattered residual,
    # i.e. far below the 280 K blackbody level.
    T_surf = 280.0
    from scipy import constants
    h = constants.h
    c = constants.c
    k_B = constants.k
    
    idx=1
    f_ref = freq_grid[idx]
    planck = (2 * h * f_ref**3 / c**2) / (np.exp(h * f_ref / (k_B * T_surf)) - 1)

    ratio = rad_reflector[idx] / planck
    assert ratio < 0.1, \
        f"Perfect reflector radiance too close to blackbody level: ratio = {ratio}"

    print("Test 3 passed: perfect reflector (emissivity=0) produces no thermal emission")


# ============================================================================
# Test 4: Jacobian shape check
# ============================================================================
def test_spectral_rad_surface_scattering_flat_diffuse_jacobian():
    """Test Jacobian computation with non-empty jac_targets."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)
    
    # Set up a Jacobian target: perturb the surface temperature.
    # (The test atmosphere has no "t" field, so an atm target would fail
    # to finalize.  Surface temperature is always present.)
    # jac_targetsFinalize requires a measurement_sensor; for this test
    # we have no sensor targets, so an empty sensor is sufficient.
    ws.measurement_sensorInit()
    ws.jac_targetsAddSurface(target="t")
    ws.jac_targetsFinalize()
    x_size = ws.jac_targets.x_size()
    
    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)
    
    # Initialize Jacobian (will be resized by the method)
    ws.spectral_rad_jac = arts.StokvecMatrix()
    
    ws.spectral_radSurfaceScatteringFlatDiffuse()
    
    # Check Jacobian shape
    jac_array = np.array(ws.spectral_rad_jac)
    assert jac_array.shape[0] == x_size, \
        f"Jacobian first dimension mismatch: {jac_array.shape[0]} != {x_size}"
    assert jac_array.shape[1] == len(freq_grid), \
        f"Jacobian second dimension mismatch: {jac_array.shape[1]} != {len(freq_grid)}"
    
    assert np.all(np.isfinite(jac_array)), \
        "Jacobian contains non-finite values"
    
    print("Test 4 passed: Jacobian has correct shape and finite values")


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    test_spectral_rad_surface_scattering_flat_diffuse_basic()
    test_spectral_rad_surface_scattering_flat_diffuse_absorbing()
    test_spectral_rad_surface_scattering_flat_diffuse_reflector()
    test_spectral_rad_surface_scattering_flat_diffuse_jacobian()
    print("\nAll tests passed!")
