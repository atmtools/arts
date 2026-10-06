"""Tests for spectral_radSurfaceScatteringInit and chained surface scattering agendas.

The *spectral_radSurfaceScattering* methods add their contribution to
*spectral_rad* and *spectral_rad_jac*, which are sized and zeroed by
*spectral_radSurfaceScatteringInit*.  Chaining several of them inside one
agenda must therefore give the sum of the individual contributions.

Verifies:
1. spectral_radSurfaceScatteringInit sizes and zeroes both outputs
2. A scattering method without the init call raises
3. Chaining Diffuse + Specular adds radiance and jacobian
4. Chaining Diffuse + DiffuseDirect adds radiance and jacobian
5. Chaining a specular method onto a pure diffuse model leaves the diffuse
   result untouched
"""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts

SUN_DISTANCE = 1.496e11
SUN_RADIUS = 6.957e8


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
    """Create and configure a minimal workspace for chained surface scattering."""
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

    # Ray point at the surface, looking straight down at the surface
    ws.ray_point = arts.PropagationPathPoint()
    ws.ray_point.pos = [0.0, 0.0, 0.0]
    ws.ray_point.los = [180.0, 0.0]

    # Single unit-weight quadrature direction, straight up and above the horizon
    ws.zen_grid = arts.ZenGrid([0.0])
    ws.az_grid = arts.AziGrid([0.0])
    ws.zen_grid_weights = arts.Vector([1.0])
    ws.az_grid_weights = arts.Vector([1.0])

    if suns is not None:
        ws.suns = suns

    # Agendas
    set_cmb_incoming_agenda(ws)
    ws.ray_path_observer_agendaSetGeometric()

    return ws


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


def flat_scalar_models(freq_grid, reflectivity):
    """MapOfSurfaceScatteringModel with a specular-only flat-scalar model."""
    field = arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[freq_grid],
        data=np.full(len(freq_grid), reflectivity).tolist(),
    )
    models = arts.MapOfSurfaceScatteringModel()
    models.add("flat_scalar", arts.FlatScalarSurfaceScatterer(field))
    return models


def blended_models(freq_grid, r_diffuse, r_specular):
    """Diffuse and specular models present at the same surface point."""
    models = lambertian_models(freq_grid, r_diffuse)
    field = arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[freq_grid],
        data=np.full(len(freq_grid), r_specular).tolist(),
    )
    models.add("flat_scalar", arts.FlatScalarSurfaceScatterer(field))
    return models


def add_surface_mask(ws, tag_key):
    """Add the surface property tag mask to ws.surf_field."""
    ws.surf_field[arts.SurfacePropertyTag(tag_key)] = 1.0


def stokes_array(ws, nf):
    return np.array([[float(ws.spectral_rad[i][s]) for s in range(4)] for i in range(nf)])


def jac_array(ws):
    return np.array(ws.spectral_rad_jac)


# ============================================================================
# Test 1: The init method sizes and zeroes the outputs
# ============================================================================
def test_init_sizes_and_zeroes():
    """spectral_radSurfaceScatteringInit must size and zero both outputs."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    ws.spectral_rad = arts.StokvecVector([arts.Stokvec([1.0, 1.0, 1.0, 1.0]) for _ in freq_grid])
    ws.spectral_radSurfaceScatteringInit()

    assert len(ws.spectral_rad) == len(freq_grid), \
        f"spectral_rad size mismatch: {len(ws.spectral_rad)} != {len(freq_grid)}"

    rad = stokes_array(ws, len(freq_grid))
    assert np.all(rad == 0.0), f"spectral_rad not zeroed by init: {rad}"

    jac = jac_array(ws)
    assert jac.shape[0] == 0 and jac.shape[1] == len(freq_grid), \
        f"spectral_rad_jac shape mismatch for empty jac_targets: {jac.shape}"

    print("Test 1 passed: init sizes and zeroes spectral_rad and spectral_rad_jac")


# ============================================================================
# Test 2: The scattering methods require the init call
# ============================================================================
def test_methods_require_init():
    """Calling a scattering method without the init method must raise."""
    freq_grid = [10e9, 100e9, 183e9]

    for name in ["Diffuse", "Specular", "DiffuseDirect", "SpecularDirect"]:
        ws = setup_workspace_base(freq_grid, suns=[make_sun(0.0, 0.0)])
        add_surface_mask(ws, "lambertian")
        ws.surface_models = lambertian_models(freq_grid, 0.5)

        try:
            getattr(ws, f"spectral_radSurfaceScattering{name}")()
        except RuntimeError as error:
            assert "spectral_radSurfaceScatteringInit" in str(error), \
                f"{name}: error does not name the init method:\n{error}"
        else:
            raise AssertionError(f"{name}: expected RuntimeError without the init call")

    print("Test 2 passed: scattering methods demand spectral_radSurfaceScatteringInit")


# ============================================================================
# Test 3: Chaining Diffuse + Specular adds both contributions
# ============================================================================
def test_chain_diffuse_plus_specular():
    """A diffuse and a specular model at the same point must add up when the
    two methods are chained after a single init call."""
    freq_grid = [10e9, 100e9, 183e9]

    def build():
        ws = setup_workspace_base(freq_grid)
        add_surface_mask(ws, "lambertian")
        add_surface_mask(ws, "flat_scalar")
        ws.surface_models = blended_models(freq_grid, 0.5, 0.3)
        return ws

    ws = build()
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()
    rad_diffuse = stokes_array(ws, len(freq_grid))

    ws = build()
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringSpecular()
    rad_specular = stokes_array(ws, len(freq_grid))

    ws = build()
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()
    ws.spectral_radSurfaceScatteringSpecular()
    rad_chain = stokes_array(ws, len(freq_grid))

    expected = rad_diffuse + rad_specular
    assert np.allclose(rad_chain, expected, rtol=1e-12, atol=0.0), \
        f"Chained agenda did not add up:\n  chain    = {rad_chain}\n  expected = {expected}"

    print("Test 3 passed: chained diffuse + specular equals the sum of both")


# ============================================================================
# Test 4: Chaining adds the jacobians as well
# ============================================================================
def test_chain_diffuse_plus_specular_jacobian():
    """The chained jacobian must equal the sum of the individual jacobians."""
    freq_grid = [10e9, 100e9, 183e9]

    def build():
        ws = setup_workspace_base(freq_grid)
        ws.measurement_sensorInit()
        ws.jac_targetsAddSurface(target="t")
        ws.jac_targetsFinalize()
        add_surface_mask(ws, "lambertian")
        add_surface_mask(ws, "flat_scalar")
        ws.surface_models = blended_models(freq_grid, 0.5, 0.3)
        return ws

    ws = build()
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()
    jac_diffuse = jac_array(ws)

    ws = build()
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringSpecular()
    jac_specular = jac_array(ws)

    ws = build()
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()
    ws.spectral_radSurfaceScatteringSpecular()
    jac_chain = jac_array(ws)

    assert jac_chain.shape == jac_diffuse.shape == jac_specular.shape, \
        f"Jacobian shape mismatch: {jac_chain.shape}, {jac_diffuse.shape}, {jac_specular.shape}"

    expected = jac_diffuse + jac_specular
    assert np.allclose(jac_chain, expected, rtol=1e-10, atol=0.0), \
        f"Chained jacobian did not add up:\n  chain    = {jac_chain}\n  expected = {expected}"

    print("Test 4 passed: chained jacobian equals the sum of both jacobians")


# ============================================================================
# Test 5: Chaining Diffuse + DiffuseDirect adds both contributions
# ============================================================================
def test_chain_diffuse_plus_diffuse_direct():
    """A quadrature-based and a sun-beam-based diffuse method must add up."""
    freq_grid = [10e9, 100e9, 183e9]
    suns = [make_sun(0.0, 0.0)]

    def build():
        ws = setup_workspace_base(freq_grid, suns=suns)
        add_surface_mask(ws, "lambertian")
        ws.surface_models = lambertian_models(freq_grid, 0.5)
        return ws

    ws = build()
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()
    rad_diffuse = stokes_array(ws, len(freq_grid))

    ws = build()
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuseDirect()
    rad_direct = stokes_array(ws, len(freq_grid))

    ws = build()
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()
    ws.spectral_radSurfaceScatteringDiffuseDirect()
    rad_chain = stokes_array(ws, len(freq_grid))

    expected = rad_diffuse + rad_direct
    assert np.allclose(rad_chain, expected, rtol=1e-12, atol=0.0), \
        f"Chained agenda did not add up:\n  chain    = {rad_chain}\n  expected = {expected}"

    print("Test 5 passed: chained diffuse + diffuse-direct equals the sum of both")


# ============================================================================
# Test 6: A specular method on a pure diffuse model adds nothing
# ============================================================================
def test_chain_specular_on_diffuse_only_model():
    """Specular methods contribute nothing for a diffuse-only model, so the
    chained result must be identical to the diffuse result."""
    freq_grid = [10e9, 100e9, 183e9]
    suns = [make_sun(0.0, 0.0)]

    def build():
        ws = setup_workspace_base(freq_grid, suns=suns)
        add_surface_mask(ws, "lambertian")
        ws.surface_models = lambertian_models(freq_grid, 0.5)
        return ws

    ws = build()
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()
    rad_diffuse = stokes_array(ws, len(freq_grid))

    ws = build()
    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()
    ws.spectral_radSurfaceScatteringSpecular()
    ws.spectral_radSurfaceScatteringSpecularDirect()
    rad_chain = stokes_array(ws, len(freq_grid))

    assert np.allclose(rad_chain, rad_diffuse, rtol=1e-12, atol=0.0), \
        f"Specular methods changed a diffuse-only result:\n  chain   = {rad_chain}\n  diffuse = {rad_diffuse}"

    print("Test 6 passed: specular methods add nothing to a diffuse-only model")


if __name__ == "__main__":
    test_init_sizes_and_zeroes()
    test_methods_require_init()
    test_chain_diffuse_plus_specular()
    test_chain_diffuse_plus_specular_jacobian()
    test_chain_diffuse_plus_diffuse_direct()
    test_chain_specular_on_diffuse_only_model()
