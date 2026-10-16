"""Tests for the sun-beam exclusion gate (``exclude_suns``) of the non-Direct
*spectral_radSurfaceScattering* methods.

Without the gate, the predefined ``SurfaceScatteringModel`` option counts the sun
twice: *spectral_radSurfaceScatteringDiffuse* / *…Specular* pick the sun up through
*spectral_rad_incoming_agenda* along a quadrature / mirror direction that contains
the solar disc, while *…DiffuseDirect* / *…SpecularDirect* add the same sun again
as a delta beam.  With ``exclude_suns = 1`` such directions are dropped entirely
(no trace, no contribution), so the ``Direct`` methods count each sun exactly once.

The gate is purely geometric (the same ``hit_sun`` disc test as the first gate of
*…SpecularDirect*), so the incoming agenda is the constant cosmic-microwave-background
one and the tests are about the geometry, not about sun radiance.

Verifies:
1. Specular gate, closed form: exclude_suns = 0 keeps R*I_CMB + eps*B(T_surf),
   exclude_suns = 1 drops the reflected term and keeps the emission term
2. Specular gate with no emission: R = 1 (eps = 0) and exclude_suns = 1 give
   exactly zero radiance
3. Specular de-duplication in the full chain: with R = 1 the full option equals
   the DirectOnly option exactly (before the fix it was exactly twice)
4. Diffuse gate, closed form: single unit-weight up-looking quadrature direction
   that is exactly the sun; exclude_suns gates it
5. Diffuse de-duplication in the full chain: with r = 1 the full option equals
   the DirectOnly option exactly
6. No regression when nothing is hit: sun at 45 deg, exclude_suns = 0 and 1 give
   identical results, and full == DiffuseOnly + DirectOnly - emission holds
7. Jacobian sanity: with a surface target and eps = 0 models the jacobian keeps
   its shape and is exactly zero

The ray_point.los uses the *upward* propagation convention: a nadir path stores
los = [0, 180] at the surface point, so the flat-surface mirror direction is
[0, 0] straight up -- exactly where make_sun(0, 0) sits.
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
    """Create and configure a minimal workspace for the sun-gate tests."""
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

    # Ray point at the surface with the upward propagation convention: a nadir
    # path stores los = [0, 180] at the surface point, so the flat-surface
    # mirror direction is the looking direction [0, 0] straight up.
    ws.ray_point = arts.PropagationPathPoint()
    ws.ray_point.pos = [0.0, 0.0, 0.0]
    ws.ray_point.los = [0.0, 180.0]

    # Single unit-weight quadrature direction, straight up -- the same direction
    # as the specular mirror direction and the sun for the default sun position
    ws.zen_grid = arts.ZenGrid([0.0])
    ws.az_grid = arts.AziGrid([0.0])
    ws.zen_grid_weights = arts.Vector([1.0])
    ws.az_grid_weights = arts.Vector([1.0])

    # Suns: default is a single sun centred exactly at the zenith
    if suns is None:
        suns = [make_sun(0.0, 0.0)]
    ws.suns = suns

    # Agendas -- CMB incoming gives a known constant radiance independent of
    # the sun; the geometric observer agenda makes the refractive LOS search
    # of the Direct methods reduce to the geometric LOS
    set_cmb_incoming_agenda(ws)
    ws.ray_path_observer_agendaSetGeometric()

    return ws


def flat_scalar_models(freq_grid, reflectivity):
    """MapOfSurfaceScatteringModel with a specular-only FlatScalar model."""
    field = arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[freq_grid],
        data=np.full(len(freq_grid), reflectivity).tolist(),
    )
    models = arts.MapOfSurfaceScatteringModel()
    models.add("flat_scalar", arts.FlatScalarSurfaceScatterer(field))
    return models


def lambertian_models(freq_grid, reflectivity):
    """MapOfSurfaceScatteringModel with a diffuse-only Lambertian model."""
    field = arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[freq_grid],
        data=np.full(len(freq_grid), reflectivity).tolist(),
    )
    models = arts.MapOfSurfaceScatteringModel()
    models.add("lambertian", arts.LambertianSurfaceScatterer(field))
    return models


def add_surface_mask(ws, tag_key):
    """Add the surface property tag mask to ws.surf_field."""
    ws.surf_field[arts.SurfacePropertyTag(tag_key)] = 1.0


def stokes_array(ws, nf):
    return np.array([[float(ws.spectral_rad[i][s]) for s in range(4)] for i in range(nf)])


def jac_array(ws):
    return np.array(ws.spectral_rad_jac)


def run_method(freq_grid, method, models, tag_key, exclude_suns, suns=None):
    """Init + a single non-Direct method with the given exclude_suns gin."""
    ws = setup_workspace_base(freq_grid, suns=suns)
    add_surface_mask(ws, tag_key)
    ws.surface_models = models
    ws.spectral_radSurfaceScatteringInit()
    getattr(ws, f"spectral_radSurfaceScattering{method}")(exclude_suns=exclude_suns)
    return stokes_array(ws, len(freq_grid))


def run_option(freq_grid, option, models, tag_key, suns=None):
    """Run a predefined spectral_rad_surface_agenda option, return (rad, jac)."""
    ws = setup_workspace_base(freq_grid, suns=suns)
    add_surface_mask(ws, tag_key)
    ws.surface_models = models
    ws.spectral_rad_surface_agendaSet(option=option)
    ws.spectral_rad_surface_agendaExecute()
    return stokes_array(ws, len(freq_grid)), jac_array(ws)


def cmb(nf):
    return np.array([planck(f, T_CMB) for f in [10e9, 100e9, 183e9]])[:nf]


# ============================================================================
# Test 1: Specular gate, closed form
# ============================================================================
def test_specular_gate_closed_form():
    """FlatScalar R = 0.5 with the sun exactly in the mirror direction.

    exclude_suns = 0 traces the mirror direction and reflects the CMB radiance
    there (today's behaviour); exclude_suns = 1 drops the direction entirely,
    leaving exactly the emission term.
    """
    freq_grid = [10e9, 100e9, 183e9]
    nf = len(freq_grid)
    R = 0.5

    rad_in = run_method(freq_grid, "Specular", flat_scalar_models(freq_grid, R), "flat_scalar", 0)
    expected_in = np.array([R * planck(f, T_CMB) + (1.0 - R) * planck(f, T_SURF) for f in freq_grid])
    assert np.allclose(rad_in[:, 0], expected_in, rtol=1e-7, atol=0.0), (
        f"exclude_suns=0 must keep the mirror-direction radiation:\n  got      = {rad_in[:, 0]}\n  expected = {expected_in}"
    )

    rad_out = run_method(freq_grid, "Specular", flat_scalar_models(freq_grid, R), "flat_scalar", 1)
    expected_out = np.array([(1.0 - R) * planck(f, T_SURF) for f in freq_grid])
    assert np.allclose(rad_out[:, 0], expected_out, rtol=1e-7, atol=0.0), (
        f"exclude_suns=1 must drop the reflected term and keep the emission:\n  got      = {rad_out[:, 0]}\n  expected = {expected_out}"
    )

    assert np.all(rad_out[:, 1:] == 0.0), "Stokes Q/U/V must stay zero"

    print("Test 1 passed: specular gate switches the reflected term on and off")


# ============================================================================
# Test 2: Specular gate with no emission
# ============================================================================
def test_specular_gate_no_emission():
    """R = 1 (eps = 0): exclude_suns = 1 must give exactly zero radiance.

    This isolates the sun double count from the emission term (which is
    identically zero here).
    """
    freq_grid = [10e9, 100e9, 183e9]

    rad = run_method(freq_grid, "Specular", flat_scalar_models(freq_grid, 1.0), "flat_scalar", 1)
    assert np.all(rad == 0.0), f"Gated pure reflector must be exactly zero:\n{rad}"

    print("Test 2 passed: gated pure specular reflector is exactly zero")


# ============================================================================
# Test 3: Specular de-duplication in the full chain
# ============================================================================
def test_specular_full_chain_dedup():
    """R = 1, sun at zenith: the full option must equal DirectOnly exactly.

    The full option runs Specular with exclude_suns = 1, so the zenith mirror
    direction is dropped and the sun enters only through SpecularDirect.
    Before the fix this assertion measured exactly 2 * I_CMB.
    """
    freq_grid = [10e9, 100e9, 183e9]
    models = flat_scalar_models(freq_grid, 1.0)

    rad_full, _ = run_option(freq_grid, "SurfaceScatteringModel", models, "flat_scalar")
    rad_direct, _ = run_option(freq_grid, "SurfaceScatteringModelDirectOnly", models, "flat_scalar")

    expected = np.array([planck(f, T_CMB) for f in freq_grid])
    assert np.allclose(rad_direct[:, 0], expected, rtol=1e-7, atol=0.0), (
        f"DirectOnly baseline must be the CMB beam:\n  got      = {rad_direct[:, 0]}\n  expected = {expected}"
    )
    assert np.allclose(rad_full, rad_direct, rtol=1e-12, atol=0.0), (
        f"Full option must count the sun exactly once:\n  full     = {rad_full[:, 0]}\n  direct   = {rad_direct[:, 0]}"
    )

    print("Test 3 passed: full chain equals DirectOnly for a pure specular reflector")


# ============================================================================
# Test 4: Diffuse gate, closed form
# ============================================================================
def test_diffuse_gate_closed_form():
    """Lambertian r = 0.5, single unit-weight up-looking quadrature direction.

    The only quadrature direction is exactly the sun, so exclude_suns = 0 gives
    r*I_CMB + eps*B(T_surf) and exclude_suns = 1 gives eps*B(T_surf).
    """
    freq_grid = [10e9, 100e9, 183e9]
    r = 0.5

    rad_in = run_method(freq_grid, "Diffuse", lambertian_models(freq_grid, r), "lambertian", 0)
    expected_in = np.array([r * planck(f, T_CMB) + (1.0 - r) * planck(f, T_SURF) for f in freq_grid])
    assert np.allclose(rad_in[:, 0], expected_in, rtol=1e-7, atol=0.0), (
        f"exclude_suns=0 must keep the sun-containing quadrature cell:\n  got      = {rad_in[:, 0]}\n  expected = {expected_in}"
    )

    rad_out = run_method(freq_grid, "Diffuse", lambertian_models(freq_grid, r), "lambertian", 1)
    expected_out = np.array([(1.0 - r) * planck(f, T_SURF) for f in freq_grid])
    assert np.allclose(rad_out[:, 0], expected_out, rtol=1e-7, atol=0.0), (
        f"exclude_suns=1 must drop the sun cell and keep the emission:\n  got      = {rad_out[:, 0]}\n  expected = {expected_out}"
    )

    print("Test 4 passed: diffuse gate switches the sun cell on and off")


# ============================================================================
# Test 5: Diffuse de-duplication in the full chain
# ============================================================================
def test_diffuse_full_chain_dedup():
    """r = 1 (eps = 0), sun at zenith: full option == DirectOnly exactly.

    The single quadrature direction is the sun; with exclude_suns = 1 the
    Diffuse method drops it and the sun enters only through DiffuseDirect.
    Before the fix this assertion measured exactly 2 * I_CMB.
    """
    freq_grid = [10e9, 100e9, 183e9]
    models = lambertian_models(freq_grid, 1.0)

    rad_full, _ = run_option(freq_grid, "SurfaceScatteringModel", models, "lambertian")
    rad_direct, _ = run_option(freq_grid, "SurfaceScatteringModelDirectOnly", models, "lambertian")

    expected = np.array([planck(f, T_CMB) for f in freq_grid])
    assert np.allclose(rad_direct[:, 0], expected, rtol=1e-7, atol=0.0), (
        f"DirectOnly baseline must be the CMB beam:\n  got      = {rad_direct[:, 0]}\n  expected = {expected}"
    )
    assert np.allclose(rad_full, rad_direct, rtol=1e-12, atol=0.0), (
        f"Full option must count the sun exactly once:\n  full     = {rad_full[:, 0]}\n  direct   = {rad_direct[:, 0]}"
    )

    print("Test 5 passed: full chain equals DirectOnly for a pure diffuse reflector")


# ============================================================================
# Test 6: No regression when nothing is hit
# ============================================================================
def test_no_gate_effect_when_sun_not_hit():
    """Sun at 45 deg: the gate must not change anything.

    Neither the up-looking quadrature direction nor the zenith mirror direction
    contains the sun, so exclude_suns = 0 and 1 agree.  The full option runs the
    Direct methods with include_emission = 0, so it carries the sub-surface
    emission once while DiffuseOnly and DirectOnly each carry it once -- the
    additivity identity full == DiffuseOnly + DirectOnly - emission holds for
    both model types.
    """
    freq_grid = [10e9, 100e9, 183e9]
    suns = [make_sun(45.0, 0.0)]
    r = 0.5

    for tag_key, models in [("flat_scalar", flat_scalar_models(freq_grid, r)),
                            ("lambertian", lambertian_models(freq_grid, r))]:
        method = "Specular" if tag_key == "flat_scalar" else "Diffuse"
        rad_off = run_method(freq_grid, method, models, tag_key, 0, suns=suns)
        rad_on = run_method(freq_grid, method, models, tag_key, 1, suns=suns)
        assert np.allclose(rad_off, rad_on, rtol=1e-12, atol=0.0), (
            f"{tag_key}: gate must be inert when no direction hits the sun:\n  off = {rad_off[:, 0]}\n  on  = {rad_on[:, 0]}"
        )

        rad_full, jac_full = run_option(freq_grid, "SurfaceScatteringModel", models, tag_key, suns=suns)
        rad_diffuse, jac_diffuse = run_option(freq_grid, "SurfaceScatteringModelDiffuseOnly", models, tag_key, suns=suns)
        rad_direct, jac_direct = run_option(freq_grid, "SurfaceScatteringModelDirectOnly", models, tag_key, suns=suns)
        # DirectOnly without suns is pure sub-surface emission -- the exact term
        # the full option removes from the Direct methods via include_emission = 0
        rad_emis, _ = run_option(freq_grid, "SurfaceScatteringModelDirectOnly", models, tag_key, suns=[])

        expected = rad_diffuse + rad_direct - rad_emis
        assert np.allclose(rad_full, expected, rtol=1e-10, atol=0.0), (
            f"{tag_key}: full != DiffuseOnly + DirectOnly - emission with the sun off-axis:\n  full     = {rad_full[:, 0]}\n  expected = {expected[:, 0]}"
        )

    print("Test 6 passed: gate inert and emission-corrected additivity intact when no direction hits the sun")


# ============================================================================
# Test 7: Jacobian sanity
# ============================================================================
def test_jacobian_zero_with_pure_reflectors():
    """With eps = 0 models and the CMB background the jacobian is exactly zero.

    jac_targetsAddSurface(target='t') gives a non-empty jacobian of shape
    (x_size, nf, 4); the emission term vanishes (eps = 0) and the CMB incoming
    agenda carries no jacobian.
    """
    freq_grid = [10e9, 100e9, 183e9]
    nf = len(freq_grid)

    for tag_key, models in [("flat_scalar", flat_scalar_models(freq_grid, 1.0)),
                            ("lambertian", lambertian_models(freq_grid, 1.0))]:
        ws = setup_workspace_base(freq_grid)
        ws.measurement_sensorInit()
        ws.jac_targetsAddSurface(target="t")
        ws.jac_targetsFinalize()
        x_size = ws.jac_targets.x_size()
        assert x_size > 0, "Expected a non-empty jacobian with a surface target"

        add_surface_mask(ws, tag_key)
        ws.surface_models = models
        ws.spectral_radSurfaceScatteringInit()
        method = "Specular" if tag_key == "flat_scalar" else "Diffuse"
        getattr(ws, f"spectral_radSurfaceScattering{method}")(exclude_suns=1)

        jac = jac_array(ws)
        assert jac.shape == (x_size, nf, 4), f"{tag_key}: jacobian shape mismatch: {jac.shape} != {(x_size, nf, 4)}"
        assert np.all(jac == 0.0), f"{tag_key}: jacobian must be exactly zero:\n{jac}"

    print("Test 7 passed: jacobian keeps its shape and is exactly zero for eps = 0")


if __name__ == "__main__":
    test_specular_gate_closed_form()
    test_specular_gate_no_emission()
    test_specular_full_chain_dedup()
    test_diffuse_gate_closed_form()
    test_diffuse_full_chain_dedup()
    test_no_gate_effect_when_sun_not_hit()
    test_jacobian_zero_with_pure_reflectors()
    print("\nAll tests passed!")
