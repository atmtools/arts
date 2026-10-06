"""Tests for spectral_radSurfaceScatteringSpecularDirect workspace method.

Verifies:
1. Smoke test: finite, non-negative output
2. Pure diffuse (Lambertian) model: specular tensors empty -> zero output
3. Sun centred in the glint direction: closed form R*I_CMB + (1-R)*B(T_surf)
4. Sun-extent rule: sun just outside the solar disc -> zero contribution; same
   offset with sun.radius inflated so beta <= alpha -> non-zero contribution
5. Sun below the local surface horizon -> emission only
6. No suns -> emission only
7. Multi-sun: one in the glint direction + one outside -> single contributing sun
8. Tilted surface: glint appears only for the sun mirrored by the actual
   surface normal into the ray direction
9. FresnelSurfaceScatterer at normal incidence: closed form with
   R = ((n1 - n2)/(n1 + n2))^2
10. Jacobian shape with jac_targetsAddSurface(target="t")

The ray_point.los uses the *upward* propagation convention: a nadir path
stores los = [0, 180] at the surface point, so specular_losNormal gives the
looking direction [0, 0] straight up for a flat surface.
"""

import numpy as np
import pyarts3 as pyarts
from scipy import constants

arts = pyarts.arts

T_CMB = 2.725  # Constant::cosmic_microwave_background_temperature

SUN_DISTANCE = 1.496e11
SUN_RADIUS = 6.957e8

T_SURF = 280.0


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
    """Incoming agenda producing the uniform cosmic microwave background."""
    @pyarts.workspace.arts_agenda(ws=ws, fix=True)
    def spectral_rad_incoming_agenda(ws):
        ws.spectral_radUniformCosmicBackground()
        ws.spectral_rad_jacEmpty()


def setup_workspace_base(freq_grid, suns=None):
    """Create and configure a minimal workspace for the specular-direct test."""
    ws = pyarts.Workspace()

    ws.freq_grid = freq_grid

    # Minimal atmosphere
    ws.abs_species = []
    ws.abs_bands = {}
    ws.spectral_propmat_agendaAuto()

    ws.atm_fieldInit(toa=100e3)

    # Surface field setup
    ws.surf_fieldEarth()
    ws.surf_field["t"] = T_SURF

    # Jacobian targets (empty for basic tests)
    ws.jac_targets = arts.JacobianTargets()

    # Ray point at the surface with the upward propagation convention:
    # a nadir path stores los = [0, 180] at the surface point, so the glint
    # direction of a flat surface is the looking direction [0, 0] straight up.
    ws.ray_point = arts.PropagationPathPoint()
    ws.ray_point.pos = [0.0, 0.0, 0.0]
    ws.ray_point.los = [0.0, 180.0]

    # Suns: default is a single sun centred exactly in the glint direction
    if suns is None:
        suns = [make_sun(0.0, 0.0)]
    ws.suns = suns

    # Agendas -- CMB incoming gives a known non-zero beam radiance; the
    # geometric observer agenda makes the refractive LOS search reduce to the
    # geometric LOS
    set_cmb_incoming_agenda(ws)
    ws.ray_path_observer_agendaSetGeometric()

    return ws


def create_surface_models(freq_grid, reflectivity):
    """Create MapOfSurfaceScatteringModel with FlatScalarSurfaceScatterer."""
    refl_data = np.full(len(freq_grid), reflectivity)

    refl_field = arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[freq_grid],
        data=refl_data.tolist()
    )

    scatterer = arts.FlatScalarSurfaceScatterer(refl_field)

    surface_models = arts.MapOfSurfaceScatteringModel()
    surface_models.add("flat_scalar", scatterer)

    return surface_models


def create_lambertian_models(freq_grid, reflectivity):
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


def create_fresnel_models(freq_grid, n2):
    """Create MapOfSurfaceScatteringModel with FresnelSurfaceScatterer (n1 = 1)."""
    n_data = np.full(len(freq_grid), n2)

    n_field = arts.SortedGriddedField1(
        name="refractive index",
        grid_names=["Frequency"],
        grids=[freq_grid],
        data=n_data.tolist()
    )

    scatterer = arts.FresnelSurfaceScatterer(n_field)

    surface_models = arts.MapOfSurfaceScatteringModel()
    surface_models.add("fresnel", scatterer)

    return surface_models


def add_surface_mask(ws, tag_key):
    """Add the surface property tag mask to ws.surf_field."""
    ws.surf_field[arts.SurfacePropertyTag(tag_key)] = 1.0


def set_tilted_surface(ws, tilt_deg=30.0):
    """Replace the surface elevation with a plane tilted in the latitudinal
    direction: h = tan(tilt) * 111320 m/deg * lat.  The surface normal at
    (0, 0) is tilted by tilt_deg towards the south, so the glint direction of
    the upward nadir ray is za = 2 * tilt_deg, aa = 180."""
    slope_m_per_deg_lat = np.tan(np.radians(tilt_deg)) * 111320.0
    lat = np.linspace(-1.0, 1.0, 5)
    lon = np.linspace(-1.0, 1.0, 5)
    h = np.repeat(slope_m_per_deg_lat * lat[:, None], len(lon), axis=1)

    ws.surf_field["h"] = arts.GeodeticField2(
        name="h",
        grid_names=["Latitude", "Longitude"],
        grids=[lat.tolist(), lon.tolist()],
        data=h.tolist(),
    )


def rad_array(ws, nf):
    return np.array([float(ws.spectral_rad[i][0]) for i in range(nf)])


def stokes_array(ws, nf):
    return np.array([[float(ws.spectral_rad[i][s]) for s in range(4)] for i in range(nf)])


def run_specular_direct(freq_grid, reflectivity, suns=None, tag_key="flat_scalar", tilt_deg=None):
    """Run the specular-direct method for a FlatScalar surface and sun list."""
    ws = setup_workspace_base(freq_grid, suns=suns)
    if tilt_deg is not None:
        set_tilted_surface(ws, tilt_deg)
    add_surface_mask(ws, tag_key)
    ws.surface_models = create_surface_models(freq_grid, reflectivity=reflectivity)
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringSpecularDirect()
    return rad_array(ws, len(freq_grid))


def emission_only(freq_grid, reflectivity):
    return np.array([(1.0 - reflectivity) * planck(f, T_SURF) for f in freq_grid])


def closed_form(freq_grid, reflectivity):
    return np.array(
        [reflectivity * planck(f, T_CMB) + (1.0 - reflectivity) * planck(f, T_SURF)
         for f in freq_grid]
    )


# ============================================================================
# Test 1: Smoke test / basic execution
# ============================================================================
def test_specular_direct_basic():
    """Test basic execution with a sun in the glint direction."""
    freq_grid = [10e9, 100e9, 183e9]
    rad = run_specular_direct(freq_grid, 0.5)

    assert len(rad) == len(freq_grid), \
        f"spectral_rad size mismatch: {len(rad)} != {len(freq_grid)}"
    assert np.all(np.isfinite(rad)), f"spectral_rad contains non-finite values: {rad}"
    assert np.all(rad >= 0), f"spectral_rad should be non-negative, got {rad}"

    print("Test 1 passed: basic execution with finite, non-negative output")


# ============================================================================
# Test 2: Pure diffuse model -> zero output
# ============================================================================
def test_specular_direct_zero_with_diffuse_model():
    """Lambertian models have zero specular BRDF and emissivity, so the
    specular-direct method must produce zero radiance."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_lambertian_models(freq_grid, reflectivity=0.5)

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringSpecularDirect()
    rad = rad_array(ws, len(freq_grid))

    assert np.all(np.isclose(rad, 0.0)), \
        f"Specular radiance should be zero for a pure diffuse model: {rad}"

    print("Test 2 passed: pure diffuse model yields zero specular radiance")


# ============================================================================
# Test 3: Closed form with the sun centred in the glint direction
# ============================================================================
def test_specular_direct_closed_form():
    """scattered == R * I_cmb, emission == (1 - R) * B(T_surf).

    The sun is centred exactly in the glint direction (straight up for the
    flat surface and the upward nadir ray).  FlatScalar: BRDF_spec = R,
    emissivity_spec = 1 - R; delta-weighted beam radiance, no cosine factor.
    """
    freq_grid = [10e9, 100e9, 183e9]
    R = 0.5

    rad = run_specular_direct(freq_grid, R)
    expected = closed_form(freq_grid, R)

    assert np.allclose(rad, expected, rtol=1e-7, atol=0.0), \
        f"Closed form violated:\n got      {rad}\n expected {expected}"

    # Pure reflector must equal the CMB exactly
    rad_1 = run_specular_direct(freq_grid, 1.0)
    expected_1 = np.array([planck(f, T_CMB) for f in freq_grid])
    assert np.allclose(rad_1, expected_1, rtol=1e-7, atol=0.0), \
        f"Pure reflector must equal CMB:\n got      {rad_1}\n expected {expected_1}"

    print("Test 3 passed: closed form R*I_CMB + (1-R)*B(T_surf)")


# ============================================================================
# Test 4: Sun-extent rule (solar disc test)
# ============================================================================
def test_specular_direct_sun_extent():
    """A sun offset by 1 deg from the glint direction is outside the solar
    disc (alpha ~ 0.2665 deg) and contributes nothing; with the sun radius
    inflated so that beta <= alpha the same geometry contributes."""
    freq_grid = [10e9, 100e9, 183e9]
    R = 0.5

    # Offset 1 deg > alpha(6.957e8) ~ 0.2665 deg -> no hit, emission only
    rad_out = run_specular_direct(freq_grid, R, suns=[make_sun(1.0, 0.0)])
    expected_out = emission_only(freq_grid, R)
    assert np.allclose(rad_out, expected_out, rtol=1e-7, atol=0.0), \
        f"Sun outside the disc must not contribute:\n got      {rad_out}\n expected {expected_out}"

    # Inflate the radius so alpha(3e9) ~ 1.146 deg > beta ~ 1 deg -> hit
    rad_in = run_specular_direct(
        freq_grid, R, suns=[make_sun(1.0, 0.0, radius=3e9)]
    )
    expected_in = closed_form(freq_grid, R)
    assert np.allclose(rad_in, expected_in, rtol=1e-7, atol=0.0), \
        f"Inflated disc must contain the glint direction:\n got      {rad_in}\n expected {expected_in}"

    print("Test 4 passed: solar disc test gates the sun contribution")


# ============================================================================
# Test 5: Sun below the local surface horizon -> emission only
# ============================================================================
def test_specular_direct_below_horizon():
    """A sun below the horizon can never satisfy beta <= alpha for the glint
    direction of an upward ray and contributes nothing."""
    freq_grid = [10e9, 100e9, 183e9]
    R = 0.5

    rad_sub = run_specular_direct(freq_grid, R, suns=[make_sun(100.0, 0.0)])
    expected = emission_only(freq_grid, R)

    assert np.allclose(rad_sub, expected, rtol=1e-7, atol=0.0), \
        f"Sub-horizon sun must not contribute:\n got      {rad_sub}\n expected {expected}"

    print("Test 5 passed: sub-horizon sun yields emission only")


# ============================================================================
# Test 6: No suns -> emission only
# ============================================================================
def test_specular_direct_no_suns():
    """An empty suns list must yield the emission term only."""
    freq_grid = [10e9, 100e9, 183e9]
    R = 0.5

    rad = run_specular_direct(freq_grid, R, suns=[])
    expected = emission_only(freq_grid, R)

    assert np.allclose(rad, expected, rtol=1e-7, atol=0.0), \
        f"No-sun case must be emission only:\n got      {rad}\n expected {expected}"

    print("Test 6 passed: no suns yields emission only")


# ============================================================================
# Test 7: Multi-sun -> only the sun in the glint direction contributes
# ============================================================================
def test_specular_direct_multi_suns():
    """One sun in the glint direction plus one outside must equal the single
    contributing sun run."""
    freq_grid = [10e9, 100e9, 183e9]
    R = 0.5

    sun_glint = make_sun(0.0, 0.0)
    sun_other = make_sun(30.0, 45.0)

    rad_one = run_specular_direct(freq_grid, R, suns=[sun_glint])
    rad_two = run_specular_direct(freq_grid, R, suns=[sun_glint, sun_other])

    assert np.allclose(rad_two, rad_one, rtol=1e-12, atol=0.0), \
        f"Non-glint sun must not contribute:\n two  {rad_two}\n one  {rad_one}"

    print("Test 7 passed: multi-sun result equals the single contributing sun")


# ============================================================================
# Test 8: Tilted surface -> glint follows the actual normal
# ============================================================================
def test_specular_direct_tilted_surface():
    """With the surface tilted 30 deg towards the south the glint direction of
    the upward nadir ray is za = 60, aa = 180: only a sun mirrored by the
    actual normal into the ray direction contributes."""
    freq_grid = [10e9, 100e9, 183e9]
    R = 0.5

    # Sun 60 deg south: mirrored by the tilted normal into the ray direction
    rad_hit = run_specular_direct(freq_grid, R, suns=[make_sun(-60.0, 0.0)], tilt_deg=30.0)
    expected_hit = closed_form(freq_grid, R)
    assert np.allclose(rad_hit, expected_hit, rtol=1e-7, atol=0.0), \
        f"Glint off tilted surface must contribute:\n got      {rad_hit}\n expected {expected_hit}"

    # Same sun on a flat surface: not in the (zenith) glint direction
    rad_flat = run_specular_direct(freq_grid, R, suns=[make_sun(-60.0, 0.0)])
    expected_flat = emission_only(freq_grid, R)
    assert np.allclose(rad_flat, expected_flat, rtol=1e-7, atol=0.0), \
        f"Sun outside flat-surface glint must not contribute:\n got      {rad_flat}\n expected {expected_flat}"

    # Overhead sun on the tilted surface: not in the tilted glint direction
    rad_tilt_none = run_specular_direct(freq_grid, R, suns=[make_sun(0.0, 0.0)], tilt_deg=30.0)
    assert np.allclose(rad_tilt_none, expected_flat, rtol=1e-7, atol=0.0), \
        f"Overhead sun off tilted glint must not contribute:\n got      {rad_tilt_none}\n expected {expected_flat}"

    print("Test 8 passed: tilted-surface glint follows the actual normal")


# ============================================================================
# Test 9: FresnelSurfaceScatterer at normal incidence
# ============================================================================
def test_specular_direct_fresnel():
    """Fresnel model with the sun in the zenith glint direction and the
    upward nadir ray: normal incidence, so BRDF = R * I4 with
    R = ((n1 - n2)/(n1 + n2))^2 and emissivity = (1 - R) * I4."""
    freq_grid = [10e9, 100e9, 183e9]
    n1, n2 = 1.0, 2.0
    R = ((n1 - n2) / (n1 + n2)) ** 2

    ws = setup_workspace_base(freq_grid)
    add_surface_mask(ws, "fresnel")
    ws.surface_models = create_fresnel_models(freq_grid, n2)

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringSpecularDirect()
    stokes = stokes_array(ws, len(freq_grid))

    expected = np.array(
        [[R * planck(f, T_CMB) + (1.0 - R) * planck(f, T_SURF), 0.0, 0.0, 0.0]
         for f in freq_grid]
    )

    assert np.allclose(stokes, expected, rtol=1e-7, atol=0.0), \
        f"Fresnel closed form violated:\n got      {stokes}\n expected {expected}"

    print("Test 9 passed: Fresnel model matches the normal-incidence closed form")


# ============================================================================
# Test 10: Jacobian shape check
# ============================================================================
def test_specular_direct_jacobian():
    """Test Jacobian computation with non-empty jac_targets."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    ws.measurement_sensorInit()
    ws.jac_targetsAddSurface(target="t")
    ws.jac_targetsFinalize()
    x_size = ws.jac_targets.x_size()

    add_surface_mask(ws, "flat_scalar")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringSpecularDirect()

    jac_array = np.array(ws.spectral_rad_jac)
    assert jac_array.shape[0] == x_size, \
        f"Jacobian first dimension mismatch: {jac_array.shape[0]} != {x_size}"
    assert jac_array.shape[1] == len(freq_grid), \
        f"Jacobian second dimension mismatch: {jac_array.shape[1]} != {len(freq_grid)}"

    assert np.all(np.isfinite(jac_array)), "Jacobian contains non-finite values"

    print("Test 10 passed: Jacobian has correct shape and finite values")


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    test_specular_direct_basic()
    test_specular_direct_zero_with_diffuse_model()
    test_specular_direct_closed_form()
    test_specular_direct_sun_extent()
    test_specular_direct_below_horizon()
    test_specular_direct_no_suns()
    test_specular_direct_multi_suns()
    test_specular_direct_tilted_surface()
    test_specular_direct_fresnel()
    test_specular_direct_jacobian()
    print("\nAll tests passed!")
