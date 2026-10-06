"""Tests for spectral_radSurfaceScatteringDiffuseDirect workspace method.

Verifies:
1. Basic execution (smoke test): method runs and produces finite output
2. Absorbing surface (r=0): output equals the surface blackbody emission
3. Beam normalization: with a known (cosmic background) incoming, the output
   equals r * I_cmb + (1 - r) * B(T_surf) in closed form (brdf = r, emiss = 1 - r)
4. Sub-horizon sun: sun below the surface-normal horizon is hard-zeroed
5. No suns: empty suns yields emission only
6. Multi-sun: two visible suns give the sum of the single-sun contributions
7. Jacobian shape: correct dimensions when jac_targets is non-empty
8. Consistency with Diffuse on a single unit-weight quadrature point
9. arts.sun.geometric_los / arts.sun.refractive_los vs sun_pathFromObserverAgenda
"""

import numpy as np
import pyarts3 as pyarts
from scipy import constants

arts = pyarts.arts

T_CMB = 2.725  # Constant::cosmic_microwave_background_temperature

SUN_DISTANCE = 1.496e11
SUN_RADIUS = 6.957e8


def planck(f, T):
    """Planck function [W m-2 sr-1 Hz-1]."""
    h   = constants.h
    c   = constants.c
    k_B = constants.k
    return (2 * h * f**3 / c**2) / (np.exp(h * f / (k_B * T)) - 1)


def make_sun(latitude, longitude, distance=SUN_DISTANCE, radius=SUN_RADIUS):
    """Create a Sun at the given sky position (geodetic lat/lon from planet center)."""
    sun = arts.Sun()
    sun.distance = distance
    sun.radius = radius
    sun.latitude = latitude
    sun.longitude = longitude
    return sun


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


def setup_workspace_base(freq_grid, suns=None):
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

    # Suns: default is a single sun above the horizon at (0, 0) on Earth
    if suns is None:
        suns = [make_sun(30.0, 45.0)]
    ws.suns = suns

    # Agendas -- CMB incoming gives a known non-zero beam radiance;
    # the geometric observer agenda makes the refractive LOS search reduce
    # to the geometric LOS
    set_cmb_incoming_agenda(ws)
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


def rad_array(ws, nf):
    return np.array([float(ws.spectral_rad[i][0]) for i in range(nf)])


def run_direct(freq_grid, reflectivity, suns=None):
    """Run the direct method for a given Lambertian reflectivity and sun list."""
    ws = setup_workspace_base(freq_grid, suns=suns)
    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=reflectivity)
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuseDirect()
    return rad_array(ws, len(freq_grid))


# ============================================================================
# Test 1: Smoke test / basic execution
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_direct_basic():
    """Test basic execution with uniform reflectivity."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuseDirect()

    assert len(ws.spectral_rad) == len(freq_grid), \
        f"spectral_rad size mismatch: {len(ws.spectral_rad)} != {len(freq_grid)}"

    rad = rad_array(ws, len(freq_grid))
    assert np.all(np.isfinite(rad)), f"spectral_rad contains non-finite values: {rad}"
    assert np.all(rad >= 0), f"spectral_rad should be non-negative, got {rad}"

    print("Test 1 passed: basic execution with finite, non-negative output")


# ============================================================================
# Test 2: Absorbing surface (reflectivity = 0) -> pure blackbody emission
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_direct_absorbing():
    """With r = 0 the output must equal the surface blackbody emission."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.0)

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuseDirect()
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
def test_spectral_rad_surface_scattering_diffuse_direct_beam_value():
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

    assert np.allclose(rad, expected, rtol=1e-7, atol=0.0), \
        f"Beam normalization violated:\n got      {rad}\n expected {expected}"

    # The scattered CMB term must be strictly positive (non-zero incoming)
    rad_1 = run_direct(freq_grid, 1.0)
    expected_1 = np.array([planck(f, T_CMB) for f in freq_grid])
    assert np.allclose(rad_1, expected_1, rtol=1e-7, atol=0.0), \
        f"Pure reflector must equal CMB:\n got      {rad_1}\n expected {expected_1}"

    print("Test 3 passed: scattered term equals r * I_cmb + (1-r) * B(T_surf)")


# ============================================================================
# Test 4: Sub-horizon sun is hard-zeroed (emission only)
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_direct_subhorizon():
    """A sun below the surface-normal horizon contributes nothing."""
    freq_grid = [10e9, 100e9, 183e9]
    r = 0.5

    # Sun at latitude 100 deg: geometric za ~ 100 > 90 at observer (0, 0),
    # i.e. below the (near-vertical) normal-based horizon
    rad_sub = run_direct(freq_grid, r, suns=[make_sun(100.0, 0.0)])

    T_surf = 280.0
    expected = np.array([(1.0 - r) * planck(f, T_surf) for f in freq_grid])

    assert np.allclose(rad_sub, expected, rtol=1e-7, atol=0.0), \
        f"Sub-horizon sun must not contribute:\n got      {rad_sub}\n expected {expected}"

    # A visible sun must give strictly more than the sub-horizon case
    rad_vis = run_direct(freq_grid, r, suns=[make_sun(30.0, 45.0)])
    assert np.all(rad_vis > rad_sub), \
        f"Visible sun must add scattered term:\n vis {rad_vis}\n sub {rad_sub}"

    print("Test 4 passed: sub-horizon sun yields emission only")


# ============================================================================
# Test 5: No suns -> emission only
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_direct_no_suns():
    """An empty suns list must yield the emission term only."""
    freq_grid = [10e9, 100e9, 183e9]
    r = 0.5

    rad = run_direct(freq_grid, r, suns=[])

    T_surf = 280.0
    expected = np.array([(1.0 - r) * planck(f, T_surf) for f in freq_grid])

    assert np.allclose(rad, expected, rtol=1e-7, atol=0.0), \
        f"No-sun case must be emission only:\n got      {rad}\n expected {expected}"

    print("Test 5 passed: no suns yields emission only")


# ============================================================================
# Test 6: Multi-sun -> sum of the single-sun contributions
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_direct_multi_suns():
    """Two visible suns: result equals the sum of the two single-sun runs
    (the emission term is counted once, so subtract the no-sun run)."""
    freq_grid = [10e9, 100e9, 183e9]
    r = 0.5

    sun_a = make_sun(30.0, 45.0)
    sun_b = make_sun(10.0, 60.0)

    rad_a = run_direct(freq_grid, r, suns=[sun_a])
    rad_b = run_direct(freq_grid, r, suns=[sun_b])
    rad_none = run_direct(freq_grid, r, suns=[])
    rad_two = run_direct(freq_grid, r, suns=[sun_a, sun_b])

    expected = rad_a + rad_b - rad_none
    assert np.allclose(rad_two, expected, rtol=1e-7, atol=0.0), \
        f"Multi-sun sum violated:\n got      {rad_two}\n expected {expected}"

    # Both suns above horizon: two scattered terms plus one emission
    T_surf = 280.0
    closed_form = np.array(
        [2.0 * r * planck(f, T_CMB) + (1.0 - r) * planck(f, T_surf) for f in freq_grid]
    )
    assert np.allclose(rad_two, closed_form, rtol=1e-7, atol=0.0), \
        f"Multi-sun closed form violated:\n got      {rad_two}\n expected {closed_form}"

    print("Test 6 passed: multi-sun result equals sum of single-sun runs")


# ============================================================================
# Test 7: Jacobian shape check
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_direct_jacobian():
    """Test Jacobian computation with non-empty jac_targets."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    ws.measurement_sensorInit()
    ws.jac_targetsAddSurface(target="t")
    ws.jac_targetsFinalize()
    x_size = ws.jac_targets.x_size()

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuseDirect()

    jac_array = np.array(ws.spectral_rad_jac)
    assert jac_array.shape[0] == x_size, \
        f"Jacobian first dimension mismatch: {jac_array.shape[0]} != {x_size}"
    assert jac_array.shape[1] == len(freq_grid), \
        f"Jacobian second dimension mismatch: {jac_array.shape[1]} != {len(freq_grid)}"

    assert np.all(np.isfinite(jac_array)), "Jacobian contains non-finite values"

    print("Test 7 passed: Jacobian has correct shape and finite values")


# ============================================================================
# Test 8: Agreement with Diffuse on a single unit-weight quadrature point
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_direct_agrees_with_diffuse():
    """The direct method must reproduce the diffuse method with a
    single-direction quadrature grid and unit weights (shared delta-weighted
    convention)."""
    freq_grid = [10e9, 100e9, 183e9]

    ws = setup_workspace_base(freq_grid)
    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    # Direct method -- take the internally estimated LOS from the sun for the
    # matching diffuse quadrature direction
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuseDirect()
    rad_direct = np.array([[float(ws.spectral_rad[i][s]) for s in range(4)]
                           for i in range(len(freq_grid))])

    za, aa = arts.sun.geometric_los(ws.suns[0], ws.ray_point.pos, ws.surf_field)

    # Diffuse method on the single direction with unit weights
    ws.zen_grid = arts.ZenGrid([za])
    ws.az_grid = arts.AziGrid([aa])
    ws.zen_grid_weights = arts.Vector([1.0])
    ws.az_grid_weights = arts.Vector([1.0])

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()
    rad_diffuse = np.array([[float(ws.spectral_rad[i][s]) for s in range(4)]
                            for i in range(len(freq_grid))])

    assert np.allclose(rad_direct, rad_diffuse, rtol=1e-12, atol=0.0), \
        f"Normalization mismatch between direct and diffuse:\n" \
        f" direct  {rad_direct}\n diffuse {rad_diffuse}"

    print("Test 8 passed: direct method agrees with single-point diffuse quadrature")


# ============================================================================
# Test 9: arts.sun LOS helpers vs sun_pathFromObserverAgenda
# ============================================================================
def test_arts_sun_los_helpers():
    """arts.sun.geometric_los / refractive_los must match the WSM sun path."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)
    sun = ws.suns[0]
    ws.sun = sun

    geo = arts.sun.geometric_los(sun, ws.ray_point.pos, ws.surf_field)
    ref = arts.sun.refractive_los(
        ws, sun, ws.ray_point.pos, ws.surf_field, ws.ray_path_observer_agenda
    )

    # With the geometric observer agenda the refractive search reduces to the
    # geometric LOS
    assert np.allclose(geo, ref, rtol=1e-12, atol=0.0), \
        f"Geometric and refractive LOS differ with geometric agenda:\n {geo}\n {ref}"

    # The WSM sun path starts at the observer with the light propagation
    # direction; mirroring it gives back the observer-pointing LOS
    ws.sun_pathFromObserverAgenda(pos=ws.ray_point.pos, angle_cut=0.0, refinement=1, just_hit=1)
    path_los = arts.path.mirror(np.array(ws.sun_path[0].los))

    assert np.allclose(geo, path_los, rtol=1e-12, atol=0.0), \
        f"arts.sun.geometric_los != mirrored sun_path front LOS:\n {geo}\n {path_los}"

    # The sun must be above the horizon at the observer for this setup
    assert 0.0 <= geo[0] < 90.0, f"Expected above-horizon sun, got za = {geo[0]}"

    print("Test 9 passed: arts.sun LOS helpers match sun_pathFromObserverAgenda")


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    test_spectral_rad_surface_scattering_diffuse_direct_basic()
    test_spectral_rad_surface_scattering_diffuse_direct_absorbing()
    test_spectral_rad_surface_scattering_diffuse_direct_beam_value()
    test_spectral_rad_surface_scattering_diffuse_direct_subhorizon()
    test_spectral_rad_surface_scattering_diffuse_direct_no_suns()
    test_spectral_rad_surface_scattering_diffuse_direct_multi_suns()
    test_spectral_rad_surface_scattering_diffuse_direct_jacobian()
    test_spectral_rad_surface_scattering_diffuse_direct_agrees_with_diffuse()
    test_arts_sun_los_helpers()
    print("\nAll tests passed!")
