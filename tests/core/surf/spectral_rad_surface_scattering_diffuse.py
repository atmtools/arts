"""Tests for spectral_radSurfaceScatteringDiffuse workspace method.

The method is the surface-normal-aware version of
spectral_radSurfaceScatteringFlatDiffuse: each quadrature direction is
checked against the horizon defined by the actual surface normal, and
directions below that horizon contribute nothing to the scattered term.

Verifies:
1. Basic execution (smoke test): method runs and produces finite output
2. Flat-surface equivalence: identical to FlatDiffuse when the surface is
   flat and the grid excludes the exact horizon
3. Tilted-surface sub-horizon gating: with a known uniform (CMB) incoming
   field and a perfect Lambertian reflector, the scattered term equals the
   flat result scaled by the fraction of quadrature weight above the
   tilted horizon
4. Jacobian shape: correct dimensions when jac_targets is non-empty
"""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts


def setup_workspace_base(freq_grid, nza=5, za_max=85.0):
    """Create and configure a minimal workspace for surface scattering tests."""
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

    # Angular grids for hemisphere integration.
    # za_max < 90 keeps the grid away from the exact horizon so that the
    # flat-surface comparison is not borderline.
    naa = 4
    za = np.linspace(0, za_max, nza)
    ws.zen_grid = arts.ZenGrid(za.tolist())
    ws.az_grid = arts.AziGrid(np.linspace(0, 360, naa, endpoint=False))

    # Quadrature weights: simple trapezoid in (za, aa)
    dza = za_max / (nza - 1) if nza > 1 else za_max
    daa = 360.0 / naa
    za_weights = np.sin(np.deg2rad(za)) * np.deg2rad(dza)
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


def set_cmb_incoming_agenda(ws):
    """Incoming agenda producing the uniform cosmic microwave background."""

    @pyarts.workspace.arts_agenda(ws=ws, fix=True)
    def spectral_rad_incoming_agenda(ws):
        ws.spectral_radUniformCosmicBackground()
        ws.spectral_rad_jacEmpty()


def create_surface_models(freq_grid, reflectivity):
    """Create MapOfSurfaceScatteringModel with LambertianSurfaceScatterer."""
    if np.isscalar(reflectivity):
        refl_data = np.full(len(freq_grid), reflectivity)
    else:
        refl_data = reflectivity

    refl_field = arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[freq_grid],
        data=refl_data.tolist(),
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


def set_tilted_surface(ws, slope_m_per_deg_lat=111320.0):
    """Replace the surface elevation with a plane tilted in the latitudinal
    direction: h = slope * lat.  With slope ~ 111320 m/deg the surface is
    tilted by ~45 degrees at (0, 0)."""
    lat = np.linspace(-1.0, 1.0, 5)
    lon = np.linspace(-1.0, 1.0, 5)
    h = np.repeat(slope_m_per_deg_lat * lat[:, None], len(lon), axis=1)

    ws.surf_field["h"] = arts.GeodeticField2(
        name="h",
        grid_names=["Latitude", "Longitude"],
        grids=[lat.tolist(), lon.tolist()],
        data=h.tolist(),
    )


def visible_weight_fraction(ws):
    """Fraction of the quadrature weight above the horizon defined by the
    actual surface normal, using the same rule as the C++ method:
    dot(ecef_los, ecef_normal) < 0."""
    ell = ws.surf_field.ellipsoid
    pos = [0.0, 0.0, 0.0]
    surf_point = ws.surf_field(0.0, 0.0)
    _, ecef_normal = arts.geodetic.geodetic_los2ecef(pos, surf_point.normal, ell)
    ecef_normal = np.array(ecef_normal)

    za = np.array(ws.zen_grid)
    aa = np.array(ws.az_grid)
    wz = np.array(ws.zen_grid_weights)
    wa = np.array(ws.az_grid_weights)

    total = 0.0
    visible = 0.0
    for i, z in enumerate(za):
        for k, a in enumerate(aa):
            _, ecef_los = arts.geodetic.geodetic_los2ecef(pos, [float(z), float(a)], ell)
            w = wz[i] * wa[k]
            total += w
            if np.dot(np.array(ecef_los), ecef_normal) < 0.0:
                visible += w

    assert visible > 0.0 and total > 0.0
    return visible / total


# ============================================================================
# Test 1: Smoke test / basic execution
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_basic():
    """Test basic execution with uniform reflectivity."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    ws.spectral_radSurfaceScatteringDiffuse()

    assert len(ws.spectral_rad) == len(freq_grid), \
        f"spectral_rad size mismatch: {len(ws.spectral_rad)} != {len(freq_grid)}"

    rad = rad_array(ws, len(freq_grid))
    assert np.all(np.isfinite(rad)), f"spectral_rad contains non-finite values: {rad}"
    assert np.all(rad > 0), f"spectral_rad should be positive (thermal scene), got {rad}"

    print("Test 1 passed: basic execution with finite, positive output")


# ============================================================================
# Test 2: Equivalence with FlatDiffuse on a flat surface
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_flat_equivalence():
    """On flat terrain with the grid away from the exact horizon every
    direction is above the surface-normal horizon, so the method must
    reproduce FlatDiffuse exactly."""
    freq_grid = [10e9, 100e9, 183e9]

    ws1 = setup_workspace_base(freq_grid)
    add_surface_mask(ws1, "lambertian")
    ws1.surface_models = create_surface_models(freq_grid, reflectivity=0.5)
    ws1.spectral_radSurfaceScatteringFlatDiffuse()
    rad_flat = rad_array(ws1, len(freq_grid))

    ws2 = setup_workspace_base(freq_grid)
    add_surface_mask(ws2, "lambertian")
    ws2.surface_models = create_surface_models(freq_grid, reflectivity=0.5)
    ws2.spectral_radSurfaceScatteringDiffuse()
    rad_diffuse = rad_array(ws2, len(freq_grid))

    assert np.allclose(rad_flat, rad_diffuse, rtol=1e-12, atol=0.0), \
        f"Flat equivalence violated:\n  flat    = {rad_flat}\n  diffuse = {rad_diffuse}"

    print("Test 2 passed: identical to FlatDiffuse on flat terrain")


# ============================================================================
# Test 3: Tilted-surface sub-horizon gating (closed form)
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_tilted_subhorizon():
    """On a ~45 deg tilted surface with a perfect Lambertian reflector
    (emissivity = 0) and a uniform CMB incoming field, the scattered term
    is r * I_cmb * sum(weights of visible directions).  The ratio to the
    flat-horizon result is therefore exactly the visible weight fraction,
    which must be strictly below 1."""
    freq_grid = [10e9, 100e9, 183e9]

    # Flat-horizon reference on the same tilted surface field
    ws1 = setup_workspace_base(freq_grid)
    set_tilted_surface(ws1)
    set_cmb_incoming_agenda(ws1)
    add_surface_mask(ws1, "lambertian")
    ws1.surface_models = create_surface_models(freq_grid, reflectivity=1.0)
    ws1.spectral_radSurfaceScatteringFlatDiffuse()
    rad_flat = rad_array(ws1, len(freq_grid))

    ws2 = setup_workspace_base(freq_grid)
    set_tilted_surface(ws2)
    set_cmb_incoming_agenda(ws2)
    add_surface_mask(ws2, "lambertian")
    ws2.surface_models = create_surface_models(freq_grid, reflectivity=1.0)
    ws2.spectral_radSurfaceScatteringDiffuse()
    rad_diffuse = rad_array(ws2, len(freq_grid))

    frac = visible_weight_fraction(ws2)
    assert 0.0 < frac < 1.0, f"Visible weight fraction not strictly inside (0, 1): {frac}"

    assert np.all(rad_flat > 0), f"Flat reference should be positive: {rad_flat}"
    assert np.all(rad_diffuse > 0), f"Tilted result should be positive: {rad_diffuse}"
    assert np.all(rad_diffuse < rad_flat), \
        f"Sub-horizon gating not effective:\n  flat    = {rad_flat}\n  diffuse = {rad_diffuse}"

    expected = rad_flat * frac
    assert np.allclose(rad_diffuse, expected, rtol=1e-9, atol=0.0), \
        f"Closed-form ratio violated:\n  got      = {rad_diffuse}\n  expected = {expected}"

    print(f"Test 3 passed: tilted-surface gating matches closed form (visible fraction {frac:.4f})")


# ============================================================================
# Test 4: Jacobian shape check
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_jacobian():
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

    ws.spectral_radSurfaceScatteringDiffuse()

    jac_array = np.array(ws.spectral_rad_jac)
    assert jac_array.shape[0] == x_size, \
        f"Jacobian first dimension mismatch: {jac_array.shape[0]} != {x_size}"
    assert jac_array.shape[1] == len(freq_grid), \
        f"Jacobian second dimension mismatch: {jac_array.shape[1]} != {len(freq_grid)}"

    assert np.all(np.isfinite(jac_array)), "Jacobian contains non-finite values"

    print("Test 4 passed: Jacobian has correct shape and finite values")


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    test_spectral_rad_surface_scattering_diffuse_basic()
    test_spectral_rad_surface_scattering_diffuse_flat_equivalence()
    test_spectral_rad_surface_scattering_diffuse_tilted_subhorizon()
    test_spectral_rad_surface_scattering_diffuse_jacobian()
    print("\nAll tests passed!")
