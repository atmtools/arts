"""Tests for spectral_radSurfaceScatteringSpecular workspace method.

Tests 1-3 register only a pure diffuse (Lambertian) model and verify the
plumbing; tests 5-6 add a Fresnel model to pin the mirror-direction physics:

1. Basic execution (smoke test): method runs and produces finite output
2. Zero-specular consistency: with a pure diffuse (Lambertian) model the
   specular BRDF and emissivity are zero, so the output radiance is zero
3. Jacobian shape: correct dimensions when jac_targets is non-empty
4. Unusable sun position (latitude outside [-90, 90]) -> user error
5. Fresnel closed form at the mirrored incoming direction (angle + Mueller)
6. exclude_suns gate hits exactly at the mirrored azimuth
"""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts


def make_sun(latitude, longitude, distance=1.496e11, radius=6.957e8):
    """Create a Sun at the given sky position (geodetic lat/lon from planet center)."""
    sun = arts.Sun()
    sun.distance = distance
    sun.radius = radius
    sun.latitude = latitude
    sun.longitude = longitude
    return sun


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
# Test 4: Unusable sun position is a user error
# ============================================================================
def test_spectral_rad_surface_scattering_specular_unusable_sun():
    """A sun that cannot be placed in the sky of the planet is rejected.

    The exclude_suns gate calls hit_sun directly, so malformed sun positions
    must be rejected by the ArrayOfSun workspace invariant before the method
    body runs -- otherwise they reach the sph2cart assertions inside hit_sun.
    """
    freq_grid = [10e9, 100e9, 183e9]
    ws = setup_workspace_base(freq_grid)
    add_surface_mask(ws, "lambertian")
    ws.surface_models = create_surface_models(freq_grid, reflectivity=0.5)

    unusable_suns = [make_sun(100.0, 0.0),
                     make_sun(0.0, 400.0),
                     make_sun(0.0, 0.0, distance=-1.0),
                     make_sun(0.0, 0.0, radius=-1.0)]

    for sun in unusable_suns:
        ws.suns = [sun]
        ws.spectral_radSurfaceScatteringInit()
        try:
            ws.spectral_radSurfaceScatteringSpecular(exclude_suns=1)
        except RuntimeError as error:
            assert "placed in the sky of the planet" in str(error), \
                f"Unexpected error for sun ({sun.latitude}, {sun.longitude}):\n{error}"
        else:
            raise AssertionError(
                f"Unusable sun ({sun.latitude}, {sun.longitude}, "
                f"{sun.distance}, {sun.radius}) was accepted")

    print("Test 4 passed: unusable sun positions raise a user error")


# ============================================================================
# Test 5: Fresnel closed form at the mirror direction
# ============================================================================
def fresnel_R(theta_deg, n1, n2):
    """Unpolarized Fresnel reflectance mean and difference at incidence angle.

    Mirrors physics_funcs::fresnel (Rv the p-form (n2 c1 - n1 c2)/(n2 c1 +
    n1 c2), Rh the s-form) and rtepack::fresnel_reflectance, returning
    (rmean, rdiff) = (0.5(|Rv|^2+|Rh|^2), 0.5(|Rv|^2-|Rh|^2)).
    """
    ti = np.radians(theta_deg)
    ci = np.cos(ti)
    ct = np.sqrt(1.0 - (n1 * np.sin(ti) / n2) ** 2)
    Rv = (n2 * ci - n1 * ct) / (n2 * ci + n1 * ct)
    Rh = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct)
    rv, rh = abs(Rv) ** 2, abs(Rh) ** 2
    return 0.5 * (rv + rh), 0.5 * (rv - rh)


def create_fresnel_models(freq_grid, n2):
    """MapOfSurfaceScatteringModel with a flat-spectrum FresnelSurfaceScatterer."""
    spectrum = arts.SortedGriddedField1(
        name="refractive_index",
        grid_names=["Frequency"],
        grids=[freq_grid],
        data=np.full(len(freq_grid), n2).tolist(),
    )
    surface_models = arts.MapOfSurfaceScatteringModel()
    surface_models.add("fresnel", arts.FresnelSurfaceScatterer(spectrum))
    return surface_models


def test_spectral_rad_surface_scattering_specular_fresnel_mirror():
    """Fresnel specular channel evaluated at the mirrored incoming direction.

    Ray looks at za=45, aa=90 from the surface; the method must mirror this
    about the actual surface normal, giving an incoming direction 45 degrees
    off the normal on the opposite azimuth.  With a uniform CMB incoming
    field the closed form is
      I: rmean(45) * I_cmb + (1 - rmean(45)) * B(T_surf)
      Q: rdiff(45) * (I_cmb - B(T_surf))
    which pins both the mirror angle and the Fresnel Mueller structure.
    """
    freq_grid = [10e9, 100e9, 183e9]
    n1, n2 = 1.0, 2.0
    za_out = 45.0

    ws = setup_workspace_base(freq_grid)
    ws.ray_point.los = [za_out, 90.0]

    @pyarts.workspace.arts_agenda(ws=ws, fix=True)
    def spectral_rad_incoming_agenda(ws):
        ws.spectral_radUniformCosmicBackground()
        ws.spectral_rad_jacEmpty()

    ws.surf_field[arts.SurfacePropertyTag("fresnel")] = 1.0
    ws.surface_models = create_fresnel_models(freq_grid, n2)

    ws.spectral_radSurfaceScatteringInit()
    ws.spectral_radSurfaceScatteringSpecular()

    rad = np.array([[float(ws.spectral_rad[i][s]) for s in range(4)]
                    for i in range(len(freq_grid))])

    T_CMB = 2.725
    T_surf = 280.0

    def planck(f, T):
        h, c, k = 6.62607015e-34, 299792458.0, 1.380649e-23
        return 2.0 * h * f**3 / c**2 / np.expm1(h * f / (k * T))

    rmean, rdiff = fresnel_R(za_out, n1, n2)
    expected = np.array(
        [[rmean * planck(f, T_CMB) + (1.0 - rmean) * planck(f, T_surf),
          rdiff * (planck(f, T_CMB) - planck(f, T_surf)),
          0.0, 0.0] for f in freq_grid]
    )

    assert np.allclose(rad, expected, rtol=1e-7, atol=0.0), \
        f"Fresnel mirror-direction closed form violated:\n got      {rad}\n expected {expected}"

    print("Test 5 passed: Fresnel specular channel matches mirror-direction closed form")


# ============================================================================
# Test 6: Sun-beam exclusion gate pins the mirror azimuth
# ============================================================================
def test_spectral_rad_surface_scattering_specular_mirror_azimuth():
    """The exclude_suns gate hits exactly at the mirrored azimuth.

    Ray los [45, 90] mirrors to sky direction [45, 270] at (0, 0), which is
    the geodetic position (lat 0, lon -45).  A sun there is gated out with
    exclude_suns=1 (emission only); a sun at (0, +45) is 90 degrees away and
    never gated.  This pins the azimuth part of the mirror direction.
    """
    freq_grid = [10e9, 100e9, 183e9]

    def run(sun, exclude):
        ws = setup_workspace_base(freq_grid)
        ws.ray_point.los = [45.0, 90.0]

        @pyarts.workspace.arts_agenda(ws=ws, fix=True)
        def spectral_rad_incoming_agenda(ws):
            ws.spectral_radUniformCosmicBackground()
            ws.spectral_rad_jacEmpty()

        ws.surf_field[arts.SurfacePropertyTag("fresnel")] = 1.0
        ws.surface_models = create_fresnel_models(freq_grid, 2.0)
        ws.suns = [sun]

        ws.spectral_radSurfaceScatteringInit()
        ws.spectral_radSurfaceScatteringSpecular(exclude_suns=exclude)
        return np.array([[float(ws.spectral_rad[i][s]) for s in range(4)]
                         for i in range(len(freq_grid))])

    rad_gated = run(make_sun(0.0, -45.0), exclude=1)
    rad_free = run(make_sun(0.0, -45.0), exclude=0)
    rad_other = run(make_sun(0.0, 45.0), exclude=1)

    # Gated sun in the mirror direction: strictly less radiance (CMB term gone)
    assert np.all(rad_gated[:, 0] < rad_free[:, 0]), \
        f"Sun at mirror direction was not gated:\n gated {rad_gated}\n free {rad_free}"
    # Sun 90 degrees off the mirror direction: gate must stay inert
    assert np.allclose(rad_other, rad_free, rtol=1e-12, atol=0.0), \
        f"Sun off the mirror direction was gated:\n other {rad_other}\n free {rad_free}"

    print("Test 6 passed: exclude_suns gate hits exactly at the mirrored azimuth")


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    test_spectral_rad_surface_scattering_specular_basic()
    test_spectral_rad_surface_scattering_specular_zero_with_diffuse_model()
    test_spectral_rad_surface_scattering_specular_jacobian()
    test_spectral_rad_surface_scattering_specular_unusable_sun()
    test_spectral_rad_surface_scattering_specular_fresnel_mirror()
    test_spectral_rad_surface_scattering_specular_mirror_azimuth()
    print("\nAll tests passed!")
