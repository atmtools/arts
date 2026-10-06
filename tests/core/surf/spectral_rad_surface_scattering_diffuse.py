"""Tests for spectral_radSurfaceScatteringDiffuse workspace method.

Each quadrature direction is checked against the horizon defined by the
actual surface normal, and directions below that horizon contribute nothing
to the scattered term.

Verifies:
1. Basic execution (smoke test): method runs and produces finite output
2. Closed-form hemisphere integration: with a perfect Lambertian reflector
   (emissivity = 0) and a uniform CMB incoming field, the radiance equals
   I_cmb times the total visible quadrature weight
3. Tilted-surface sub-horizon gating: the scattered term on a tilted surface
   equals the flat-surface result scaled by the visible weight fraction
4. Absorbing surface (r = 0): output equals the surface blackbody
5. Kirchhoff coupling: perfect reflector (r = 1) emits nothing and reflects
   no thermal radiation, so it produces far less radiance than the absorber
6. Jacobian shape: correct dimensions when jac_targets is non-empty
"""

import numpy as np
import pyarts3 as pyarts
from scipy import constants

arts = pyarts.arts

T_CMB = 2.725  # Constant::cosmic_microwave_background_temperature


def planck(f, T):
    """Planck function [W m-2 sr-1 Hz-1]."""
    h = constants.h
    c = constants.c
    k_B = constants.k
    return (2 * h * f**3 / c**2) / (np.exp(h * f / (k_B * T)) - 1)


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


def total_quadrature_weight(ws):
    return float(np.sum(np.array(ws.zen_grid_weights)) * np.sum(np.array(ws.az_grid_weights)))


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

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()

    assert len(ws.spectral_rad) == len(freq_grid), \
        f"spectral_rad size mismatch: {len(ws.spectral_rad)} != {len(freq_grid)}"

    rad = rad_array(ws, len(freq_grid))
    assert np.all(np.isfinite(rad)), f"spectral_rad contains non-finite values: {rad}"
    assert np.all(rad > 0), f"spectral_rad should be positive (thermal scene), got {rad}"

    print("Test 1 passed: basic execution with finite, positive output")


# ============================================================================
# Test 2: Closed-form hemisphere integration on flat terrain
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_closed_form():
    """With a perfect Lambertian reflector (r = 1, emissivity = 0) and a
    uniform CMB incoming field the method must return exactly
    I_cmb * sum(visible quadrature weights).  On flat terrain with the grid
    away from the exact horizon every direction is visible."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    set_cmb_incoming_agenda(ws)
    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=1.0)

    frac = visible_weight_fraction(ws)
    assert np.isclose(frac, 1.0, rtol=1e-12, atol=0.0), \
        f"All directions should be visible on flat terrain, got fraction {frac}"

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()
    rad = rad_array(ws, len(freq_grid))

    cmb = np.array([planck(f, T_CMB) for f in freq_grid])
    expected = cmb * total_quadrature_weight(ws)

    assert np.allclose(rad, expected, rtol=1e-12, atol=0.0), \
        f"Closed-form hemisphere integration violated:\n  got      = {rad}\n  expected = {expected}"

    print("Test 2 passed: hemisphere integration matches closed form on flat terrain")


# ============================================================================
# Test 3: Tilted-surface sub-horizon gating (closed form)
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_tilted_subhorizon():
    """On a ~45 deg tilted surface with a perfect Lambertian reflector
    (emissivity = 0) and a uniform CMB incoming field, the scattered term
    is r * I_cmb * sum(weights of visible directions).  The ratio to the
    flat-surface result is therefore exactly the visible weight fraction,
    which must be strictly below 1."""
    freq_grid = [10e9, 100e9, 183e9]

    # Flat-terrain reference with the identical setup
    ws1 = setup_workspace_base(freq_grid)
    set_cmb_incoming_agenda(ws1)
    add_surface_mask(ws1, "lambertian")
    ws1.surface_models = create_surface_models(freq_grid, reflectivity=1.0)
    ws1.spectral_radSurfaceScatteringInit()
    ws1.spectral_radSurfaceScatteringDiffuse()
    rad_flat = rad_array(ws1, len(freq_grid))

    ws2 = setup_workspace_base(freq_grid)
    set_tilted_surface(ws2)
    set_cmb_incoming_agenda(ws2)
    add_surface_mask(ws2, "lambertian")
    ws2.surface_models = create_surface_models(freq_grid, reflectivity=1.0)
    ws2.spectral_radSurfaceScatteringInit()
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
# Test 4: Absorbing surface (reflectivity = 0)
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_absorbing():
    """With reflectivity = 0 the scattering contribution vanishes and the
    output is the pure surface emission at the surface temperature."""
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)

    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.0)

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()
    rad_absorbing = rad_array(ws, len(freq_grid))

    assert np.all(np.isfinite(rad_absorbing)), \
        "Absorbing surface: spectral_rad contains non-finite values"
    assert np.all(rad_absorbing > 0), \
        "Absorbing surface: spectral_rad should be positive (surface emission)"

    T_surf = ws.surf_field["t"](0, 0)
    idx = 1
    ratio = rad_absorbing[idx] / planck(freq_grid[idx], T_surf)
    assert 0.9999999 < ratio < 1.0000001, \
        f"Planck ratio out of range: {ratio} (expect 1.0 for emissivity = 1)"

    print("Test 4 passed: absorbing surface produces the surface blackbody")


# ============================================================================
# Test 5: Kirchhoff coupling, perfect reflector (reflectivity = 1)
# ============================================================================
def test_spectral_rad_surface_scattering_diffuse_reflector():
    """Test the Kirchhoff coupling emissivity = 1 - reflectivity.

    The test atmosphere contains no absorption species, so there is no
    incoming thermal radiation for the surface to reflect.  Therefore a
    perfect reflector (r = 1, emissivity = 0) produces (near-)zero radiance,
    while the absorbing surface (r = 0) produces the full 280 K blackbody.
    """
    freq_grid = [10e9, 100e9, 183e9]

    ws1 = setup_workspace_base(freq_grid)
    add_surface_mask(ws1, "lambertian")
    ws1.surface_models = create_surface_models(freq_grid, reflectivity=0.0)
    ws1.spectral_radSurfaceScatteringInit()
    ws1.spectral_radSurfaceScatteringDiffuse()
    rad_absorbing = rad_array(ws1, len(freq_grid))

    ws2 = setup_workspace_base(freq_grid)
    add_surface_mask(ws2, "lambertian")
    ws2.surface_models = create_surface_models(freq_grid, reflectivity=1.0)
    ws2.spectral_radSurfaceScatteringInit()
    ws2.spectral_radSurfaceScatteringDiffuse()
    rad_reflector = rad_array(ws2, len(freq_grid))

    assert np.all(rad_reflector < rad_absorbing), \
        f"Emissivity = 1 - r violated: reflector {rad_reflector} >= absorber {rad_absorbing}"

    idx = 1
    ratio = rad_reflector[idx] / planck(freq_grid[idx], 280.0)
    assert ratio < 0.1, \
        f"Perfect reflector radiance too close to blackbody level: ratio = {ratio}"

    print("Test 5 passed: perfect reflector (emissivity = 0) produces no thermal emission")


# ============================================================================
# Test 6: Jacobian shape check
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

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringDiffuse()

    jac_array = np.array(ws.spectral_rad_jac)
    assert jac_array.shape[0] == x_size, \
        f"Jacobian first dimension mismatch: {jac_array.shape[0]} != {x_size}"
    assert jac_array.shape[1] == len(freq_grid), \
        f"Jacobian second dimension mismatch: {jac_array.shape[1]} != {len(freq_grid)}"

    assert np.all(np.isfinite(jac_array)), "Jacobian contains non-finite values"

    print("Test 6 passed: Jacobian has correct shape and finite values")


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    test_spectral_rad_surface_scattering_diffuse_basic()
    test_spectral_rad_surface_scattering_diffuse_closed_form()
    test_spectral_rad_surface_scattering_diffuse_tilted_subhorizon()
    test_spectral_rad_surface_scattering_diffuse_absorbing()
    test_spectral_rad_surface_scattering_diffuse_reflector()
    test_spectral_rad_surface_scattering_diffuse_jacobian()
    print("\nAll tests passed!")
