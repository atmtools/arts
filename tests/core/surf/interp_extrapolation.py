"""Tests for the interp_extrapolation frequency extrapolation modes.

The models interpolate their stored spectral data onto the simulation f_grid;
outside the stored frequency grid the interp_extrapolation mode decides what
happens.  Verifies for each mode:

1. Linear: unlimited linear extrapolation outside the stored grid
2. Nearest: edge-value clamping outside the stored grid (the default)
3. Zero: values outside the stored grid are 0
4. None: query frequencies outside the stored grid are a user error
5. The modes apply to the spatial field variants as well
6. Fresnel: Nearest equals edge evaluation, None and Zero are user errors
7. FlatScalar: Zero gives zero reflectivity outside the stored grid
"""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts

IE = arts.InterpolationExtrapolation

# Stored spectral grid and linear ramp data
f_stored = [1e10, 1e12]
r_lo, r_hi = 0.4, 0.6

# Query grid extending beyond the stored grid on both sides
f_query = arts.Vector([1e9, 5e11, 2e12])

za = arts.Vector([0.0])
aa = arts.Vector([0.0])
surf_pt = arts.SurfacePoint()


def ramp(f):
    """The stored linear ramp evaluated at frequency f."""
    return r_lo + (r_hi - r_lo) * (f - f_stored[0]) / (f_stored[1] - f_stored[0])


def lambertian_brdf_i(sc, f_grid):
    """Diffuse BRDF Stokes-I diagonal of a Lambertian model over f_grid."""
    props = sc.get_surface_scattering_model_properties(
        surf_pt, 0.0, 0.0, f_grid, za, aa, za, aa)
    return np.array(props.brdf_matrix_diffuse)[:, 0, 0, 0, 0, 0, 0]


def make_spectrum_lambertian(mode):
    spectrum = arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[f_stored],
        data=[r_lo, r_hi],
    )
    sc = arts.LambertianSurfaceScatterer(spectrum)
    sc.interp_extrapolation = mode
    return sc


# ============================================================================
# Test 1: Linear -- unlimited linear extrapolation
# ============================================================================
def test_linear():
    sc = make_spectrum_lambertian(IE.Linear)
    got = lambertian_brdf_i(sc, f_query)
    expected = np.array([ramp(f) for f in f_query])
    assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
        f"Linear extrapolation mismatch:\n got      {got}\n expected {expected}"
    print("Test 1 passed: Linear extrapolates the ramp outside the stored grid")


# ============================================================================
# Test 2: Nearest -- edge-value clamping (and that it is the default)
# ============================================================================
def test_nearest():
    sc = make_spectrum_lambertian(IE.Nearest)
    assert sc.interp_extrapolation == IE.Nearest, "Nearest must be the default"
    got = lambertian_brdf_i(sc, f_query)
    expected = np.array([r_stored for r_stored in (r_lo, ramp(5e11), r_hi)])
    assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
        f"Nearest clamping mismatch:\n got      {got}\n expected {expected}"
    print("Test 2 passed: Nearest clamps to the stored grid edge values")


# ============================================================================
# Test 3: Zero -- 0 outside the stored grid
# ============================================================================
def test_zero():
    sc = make_spectrum_lambertian(IE.Zero)
    got = lambertian_brdf_i(sc, f_query)
    expected = np.array([0.0, ramp(5e11), 0.0])
    assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
        f"Zero extrapolation mismatch:\n got      {got}\n expected {expected}"
    print("Test 3 passed: Zero gives 0 outside the stored grid")


# ============================================================================
# Test 4: None -- user error outside the stored grid
# ============================================================================
def test_none():
    sc = make_spectrum_lambertian(IE.None_)
    try:
        lambertian_brdf_i(sc, f_query)
    except RuntimeError as e:
        assert "interp_extrapolation is None" in str(e), \
            f"Unexpected error message: {e}"
    else:
        raise AssertionError("None mode must raise for frequencies outside the grid")

    # Inside the grid the same model works fine
    inside = arts.Vector([1e10, 5e11, 1e12])
    got = lambertian_brdf_i(sc, inside)
    expected = np.array([ramp(f) for f in inside])
    assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
        f"None mode inside grid mismatch:\n got      {got}\n expected {expected}"
    print("Test 4 passed: None raises a user error outside the stored grid")


# ============================================================================
# Test 5: Field variants honour the modes as well
# ============================================================================
def test_field_variants():
    field3 = arts.SortedGriddedField3(
        name="reflectivity",
        grid_names=["Latitude", "Longitude", "Frequency"],
        grids=[[-90.0, 90.0], [-180.0, 180.0], f_stored],
        data=np.full((2, 2, 2), 0.0)
             + np.array([[[r_lo, r_hi]]]),  # ramp along frequency, flat spatially
    )

    expected_nearest = np.array([r_lo, ramp(5e11), r_hi])
    expected_zero = np.array([0.0, ramp(5e11), 0.0])

    for mode, expected in [(IE.Nearest, expected_nearest), (IE.Zero, expected_zero)]:
        sc = arts.LambertianSurfaceScattererField(field3)
        sc.interp_extrapolation = mode
        got = lambertian_brdf_i(sc, f_query)
        assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
            f"Field {mode} mismatch:\n got      {got}\n expected {expected}"

    sc = arts.LambertianSurfaceScattererField(field3)
    sc.interp_extrapolation = IE.None_
    try:
        lambertian_brdf_i(sc, f_query)
    except RuntimeError as e:
        assert "interp_extrapolation is None" in str(e), \
            f"Unexpected error message: {e}"
    else:
        raise AssertionError("Field None mode must raise outside the grid")

    print("Test 5 passed: field variants honour Nearest/Zero/None")


# ============================================================================
# Test 6: Fresnel -- Nearest equals edge evaluation, None/Zero are errors
# ============================================================================
def test_fresnel_modes():
    n_lo, n_hi = 1.5, 2.5
    spectrum = arts.SortedGriddedField1(
        name="refractive_index",
        grid_names=["Frequency"],
        grids=[f_stored],
        data=[n_lo, n_hi],
    )

    def specular_i(sc, f_grid):
        props = sc.get_surface_scattering_model_properties(
            surf_pt, 0.0, 0.0, f_grid, za, aa, za, aa)
        return np.array(props.brdf_matrix_specular)[:, 0, 0, 0, 0, 0, 0]

    # Nearest at out-of-grid frequencies must equal Linear evaluated at the
    # clamped (edge) frequencies -- no Fresnel formula needed, just geometry
    sc_nearest = arts.FresnelSurfaceScatterer(spectrum)
    sc_nearest.interp_extrapolation = IE.Nearest
    got = specular_i(sc_nearest, f_query)

    sc_linear = arts.FresnelSurfaceScatterer(spectrum)
    sc_linear.interp_extrapolation = IE.Linear
    expected = specular_i(sc_linear, arts.Vector([f_stored[0], 5e11, f_stored[1]]))

    assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
        f"Fresnel Nearest != edge evaluation:\n got      {got}\n expected {expected}"

    sc_none = arts.FresnelSurfaceScatterer(spectrum)
    sc_none.interp_extrapolation = IE.None_
    try:
        specular_i(sc_none, f_query)
    except RuntimeError as e:
        assert "interp_extrapolation is None" in str(e), \
            f"Unexpected error message: {e}"
    else:
        raise AssertionError("Fresnel None mode must raise outside the grid")

    # Zero mode gives n2 = 0 outside the grid -- a user error, the refractive
    # index must stay positive
    sc_zero = arts.FresnelSurfaceScatterer(spectrum)
    sc_zero.interp_extrapolation = IE.Zero
    try:
        specular_i(sc_zero, f_query)
    except RuntimeError:
        pass
    else:
        raise AssertionError("Fresnel Zero mode must raise (n2 = 0 is invalid)")

    print("Test 6 passed: Fresnel Nearest clamps, None and Zero are user errors")


# ============================================================================
# Test 7: FlatScalar -- Zero reflectivity outside the stored grid
# ============================================================================
def test_flat_scalar_zero():
    spectrum = arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[f_stored],
        data=[r_lo, r_hi],
    )
    sc = arts.FlatScalarSurfaceScatterer(spectrum)
    sc.interp_extrapolation = IE.Zero
    props = sc.get_surface_scattering_model_properties(
        surf_pt, 0.0, 0.0, f_query, za, aa, za, aa)
    got = np.array(props.brdf_matrix_specular)[:, 0, 0, 0, 0, 0, 0]
    expected = np.array([0.0, ramp(5e11), 0.0])
    assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
        f"FlatScalar Zero mismatch:\n got      {got}\n expected {expected}"
    print("Test 7 passed: FlatScalar Zero gives zero reflectivity outside the grid")


# ============================================================================
# Main
# ============================================================================
if __name__ == "__main__":
    test_linear()
    test_nearest()
    test_zero()
    test_none()
    test_field_variants()
    test_fresnel_modes()
    test_flat_scalar_zero()
    print("\nAll interp_extrapolation tests passed!")
