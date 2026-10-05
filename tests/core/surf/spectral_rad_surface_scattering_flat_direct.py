"""Tests for spectral_radSurfaceScatteringFlatDirect workspace method.

Verifies:
1. Basic execution (smoke test): method runs and produces finite output
2. Absorbing surface (r=0): output equals the surface blackbody emission
3. Beam normalization: with a known (cosmic background) incoming, the output
   equals r * I_cmb + (1 - r) * B(T_surf) in closed form (brdf = r, emiss = 1 - r)
4. Sub-horizon beam: beam below the surface-normal horizon is hard-zeroed
5. Jacobian shape: correct dimensions when jac_targets is non-empty
6. Consistency with FlatDiffuse on a single unit-weight quadrature point
"""

import numpy as np
import pyarts3 as pyarts
from scipy import constants

arts = pyarts.arts

T_CMB = 2.725  # Constant::cosmic_microwave_background_temperature


def planck(f, T):
    """Planck function [W m-2 sr-1 Hz-1]."""
    h   = constants.h
    c   = constants.c
    k_B = constants.k
    return (2 * h * f**3 / c**2) / (np.exp(h * f / (k_B * T)) - 1)


def set_cmb_incoming_agenda(ws):
    """Incoming agenda producing the uniform cosmic microwave background.

    Uses only workspace methods (which write through the agenda output
    pointers), so it works when the agenda is invoked from C++.  The CMB
    radiance is a known non-zero incoming field, which exercises the
    scattered/BRDF path with a closed-form value.
    """
    @pyarts.workspace.arts_agenda(ws=ws, fix=True)
    def spectral_rad_incoming_agenda(ws):
        ws.spectral_radUniformCosmicBackground()
        ws.spectral_rad_jacEmpty()


def setup_workspace_base(freq_grid):
    """Create and configure a minimal workspace for direct-beam tests."""
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

    # Direct beam from above the horizon (visible at (0, 0) on Earth)
    ws.direct_beam_los = arts.Vector2([30.0, 45.0])

    # Agendas -- CMB incoming gives a known non-zero beam radiance
    set_cmb_incoming_agenda(ws)
    ws.ray_path_observer_agendaSetGeometric()

    # Initialize spectral_rad and spectral_rad_jac (will be resized by the method)
    ws.spectral_rad = arts.StokvecVector()
    ws.spectral_rad_jac = arts.StokvecMatrix()

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


def rad_array(ws, nf):
    return np.array([float(ws.spectral_rad[i][0]) for i in range(nf)])


def run_direct(freq_grid, reflectivity, direct_beam_los=(30.0, 45.0)):
    """Run the direct method for a given Lambertian reflectivity."""
    ws = setup_workspace_base(freq_grid)
    ws.direct_beam_los = arts.Vector2(list(direct_beam_los))
    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=reflectivity)
    ws.spectral_radSurfaceScatteringFlatDirect()
    return rad_array(ws, len(freq_grid))


# ============================================================================
# Test 1: Smoke test / basic execution
# ============================================================================
def test_spectral_rad_surface_scattering_flat_direct_basic():
    """Test basic execution with uniform reflectivity."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    ws.spectral_radSurfaceScatteringFlatDirect()

    assert len(ws.spectral_rad) == len(freq_grid), \
        f"spectral_rad size mismatch: {len(ws.spectral_rad)} != {len(freq_grid)}"

    rad = rad_array(ws, len(freq_grid))
    assert np.all(np.isfinite(rad)), f"spectral_rad contains non-finite values: {rad}"
    assert np.all(rad >= 0), f"spectral_rad should be non-negative, got {rad}"

    print("Test 1 passed: basic execution with finite, non-negative output")


# ============================================================================
# Test 2: Absorbing surface (reflectivity = 0) -> pure blackbody emission
# ============================================================================
def test_spectral_rad_surface_scattering_flat_direct_absorbing():
    """With r = 0 the output must equal the surface blackbody emission."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.0)

    ws.spectral_radSurfaceScatteringFlatDirect()
    rad = rad_array(ws, len(freq_grid))

    assert np.all(np.isfinite(rad)), "Absorbing surface: non-finite values"

    T_surf = ws.surf_field["t"](0, 0)
    for i, f in enumerate(freq_grid):
        ratio = rad[i] / planck(f, T_surf)
        assert 0.9999999 < ratio < 1.0000001, \
            f"Planck ratio out of range at f={f}: {ratio}"

    print("Test 2 passed: absorbing surface reproduces the surface blackbody")


# ============================================================================
# Test 3: Closed-form normalization with a known (CMB) incoming beam
# ============================================================================
def test_spectral_rad_surface_scattering_flat_direct_beam_value():
    """scattered == r * I_cmb, emission == (1 - r) * B(T_surf).

    The incoming agenda supplies the uniform cosmic microwave background, a
    known radiance.  For a Lambertian surface brdf = r and emissivity = 1 - r.
    The beam radiance is delta-weighted: no quadrature weights, no 1/pi, no
    cosine factor.  This pins the full normalization in closed form.
    """
    freq_grid = [10e9, 100e9, 183e9]
    r = 0.5

    rad = run_direct(freq_grid, r)

    T_surf = 280.0
    expected = np.array(
        [r * planck(f, T_CMB) + (1.0 - r) * planck(f, T_surf) for f in freq_grid]
    )

    assert np.allclose(rad, expected, rtol=1e-7), \
        f"Beam normalization violated:\n got      {rad}\n expected {expected}"

    # The scattered CMB term must be strictly positive (non-zero incoming)
    rad_1 = run_direct(freq_grid, 1.0)
    expected_1 = np.array([planck(f, T_CMB) for f in freq_grid])
    assert np.allclose(rad_1, expected_1, rtol=1e-7), \
        f"Pure reflector must equal CMB:\n got      {rad_1}\n expected {expected_1}"

    print("Test 3 passed: scattered term equals r * I_cmb + (1-r) * B(T_surf)")


# ============================================================================
# Test 4: Sub-horizon beam is hard-zeroed (emission only)
# ============================================================================
def test_spectral_rad_surface_scattering_flat_direct_subhorizon():
    """A beam below the surface-normal horizon contributes nothing."""
    freq_grid = [10e9, 100e9, 183e9]
    r = 0.5

    # za = 100 > 90: below the (near-vertical) normal-based horizon at (0, 0)
    rad_sub = run_direct(freq_grid, r, direct_beam_los=(100.0, 0.0))

    T_surf = 280.0
    expected = np.array([(1.0 - r) * planck(f, T_surf) for f in freq_grid])

    assert np.allclose(rad_sub, expected, rtol=1e-7), \
        f"Sub-horizon beam must not contribute:\n got      {rad_sub}\n expected {expected}"

    # A visible beam must give strictly more than the sub-horizon case
    rad_vis = run_direct(freq_grid, r, direct_beam_los=(30.0, 45.0))
    assert np.all(rad_vis > rad_sub), \
        f"Visible beam must add scattered term:\n vis {rad_vis}\n sub {rad_sub}"

    print("Test 4 passed: sub-horizon beam yields emission only")


# ============================================================================
# Test 5: Jacobian shape check
# ============================================================================
def test_spectral_rad_surface_scattering_flat_direct_jacobian():
    """Test Jacobian computation with non-empty jac_targets."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    ws.measurement_sensorInit()
    ws.jac_targetsAddSurface(target="t")
    ws.jac_targetsFinalize()
    x_size = ws.jac_targets.x_size()

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    ws.spectral_rad_jac = arts.StokvecMatrix()

    ws.spectral_radSurfaceScatteringFlatDirect()

    jac_array = np.array(ws.spectral_rad_jac)
    assert jac_array.shape[0] == x_size, \
        f"Jacobian first dimension mismatch: {jac_array.shape[0]} != {x_size}"
    assert jac_array.shape[1] == len(freq_grid), \
        f"Jacobian second dimension mismatch: {jac_array.shape[1]} != {len(freq_grid)}"

    assert np.all(np.isfinite(jac_array)), "Jacobian contains non-finite values"

    print("Test 5 passed: Jacobian has correct shape and finite values")


# ============================================================================
# Test 6: Agreement with FlatDiffuse on a single unit-weight quadrature point
# ============================================================================
def test_spectral_rad_surface_scattering_flat_direct_agrees_with_diffuse():
    """The direct method must reproduce FlatDiffuse with a single-direction
    quadrature grid and unit weights (shared delta-weighted convention)."""
    freq_grid = [10e9, 100e9, 183e9]
    za, aa = 40.0, 137.0

    ws = setup_workspace_base(freq_grid)
    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    # Direct method
    ws.direct_beam_los = arts.Vector2([za, aa])
    ws.spectral_radSurfaceScatteringFlatDirect()
    rad_direct = np.array([[float(ws.spectral_rad[i][s]) for s in range(4)]
                           for i in range(len(freq_grid))])

    # Diffuse method on the single direction with unit weights
    ws.zen_grid = arts.ZenGrid([za])
    ws.az_grid = arts.AziGrid([aa])
    ws.zen_grid_weights = arts.Vector([1.0])
    ws.az_grid_weights = arts.Vector([1.0])

    ws.spectral_rad = arts.StokvecVector()
    ws.spectral_rad_jac = arts.StokvecMatrix()
    ws.spectral_radSurfaceScatteringFlatDiffuse()
    rad_diffuse = np.array([[float(ws.spectral_rad[i][s]) for s in range(4)]
                            for i in range(len(freq_grid))])

    assert np.allclose(rad_direct, rad_diffuse, rtol=1e-12), \
        f"Normalization mismatch between direct and diffuse:\n" \
        f" direct  {rad_direct}\n diffuse {rad_diffuse}"

    print("Test 6 passed: direct method agrees with single-point FlatDiffuse")


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    test_spectral_rad_surface_scattering_flat_direct_basic()
    test_spectral_rad_surface_scattering_flat_direct_absorbing()
    test_spectral_rad_surface_scattering_flat_direct_beam_value()
    test_spectral_rad_surface_scattering_flat_direct_subhorizon()
    test_spectral_rad_surface_scattering_flat_direct_jacobian()
    test_spectral_rad_surface_scattering_flat_direct_agrees_with_diffuse()
    print("\nAll tests passed!")
