"""Tests for the per-dimension extrapolation modes of SurfaceSpectralData.

The models evaluate their SurfaceSpectralData (reflectivity or refractive
index over frequency, latitude, and longitude) on the simulation f_grid at the
surface position.  Outside the stored grids, the extrapolation mode of each
dimension and side (frq_low/frq_upp, lat_low/lat_upp, lon_low/lon_upp)
decides what happens.  Verifies for the frequency dimension:

1. Linear: unlimited linear extrapolation outside the stored grid
2. Nearest: edge-value clamping outside the stored grid
3. Zero: values outside the stored grid are 0
4. None: query frequencies outside the stored grid are a user error; this is
   the default for a multi-point frequency grid, while single-point latitude
   and longitude grids default to Nearest
5. The modes apply to spatially gridded data as well
6. Fresnel: Nearest equals edge evaluation, None and Zero are user errors
7. FlatScalar: Zero gives zero reflectivity outside the stored grid
8. The lower and upper sides of a dimension are independent
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

# Error message of a query outside the grid in None mode
NONE_ERROR = "out-of-bounds when no extrapolation is wanted"


def ramp(f):
    """The stored linear ramp evaluated at frequency f."""
    return r_lo + (r_hi - r_lo) * (f - f_stored[0]) / (f_stored[1] - f_stored[0])


def spectrum(values):
    """SurfaceSpectralData of a pure spectrum on the stored frequency grid.

    The single-point latitude and longitude grids make it constant in space.
    """
    g = arts.SpectralGeodeticField3()
    g.grids = (arts.AscendingGrid(f_stored), arts.LatGrid([0.0]), arts.LonGrid([0.0]))
    g.data = np.array(values, dtype=float).reshape(-1, 1, 1)
    return arts.SurfaceSpectralData(g)


def set_frequency_extrapolation(data, mode):
    """Set the extrapolation below and above the stored frequency grid."""
    data.frq_low = mode
    data.frq_upp = mode


def lambertian_brdf_i(sc, f_grid):
    """Diffuse BRDF Stokes-I diagonal of a Lambertian model over f_grid."""
    props = sc.get_surface_scattering_model_properties(
        surf_pt, 0.0, 0.0, f_grid, za, aa, za, aa)
    return np.array(props.brdf_matrix_diffuse)[:, 0, 0, 0, 0, 0, 0]


def make_spectrum_lambertian(mode):
    sc = arts.LambertianSurfaceScatterer(spectrum([r_lo, r_hi]))
    set_frequency_extrapolation(sc.reflectivity, mode)
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
# Test 2: Nearest -- edge-value clamping
# ============================================================================
def test_nearest():
    sc = make_spectrum_lambertian(IE.Nearest)
    assert sc.reflectivity.frq_low == IE.Nearest, "frq_low must read back Nearest"
    assert sc.reflectivity.frq_upp == IE.Nearest, "frq_upp must read back Nearest"
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
# Test 4: None -- user error outside the stored grid (and the default)
# ============================================================================
def test_none():
    # Default: None for the multi-point frequency grid, Nearest for the
    # single-point latitude and longitude grids
    sc_default = arts.LambertianSurfaceScatterer(spectrum([r_lo, r_hi]))
    refl = sc_default.reflectivity
    assert refl.frq_low == IE.None_ and refl.frq_upp == IE.None_, \
        "None must be the default for a multi-point frequency grid"
    for side in (refl.lat_low, refl.lat_upp, refl.lon_low, refl.lon_upp):
        assert side == IE.Nearest, \
            "Nearest must be the default for single-point spatial grids"

    for sc in (sc_default, make_spectrum_lambertian(IE.None_)):
        try:
            lambertian_brdf_i(sc, f_query)
        except RuntimeError as e:
            assert NONE_ERROR in str(e), f"Unexpected error message: {e}"
        else:
            raise AssertionError("None mode must raise for frequencies outside the grid")

        # Inside the grid the same model works fine
        inside = arts.Vector([1e10, 5e11, 1e12])
        got = lambertian_brdf_i(sc, inside)
        expected = np.array([ramp(f) for f in inside])
        assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
            f"None mode inside grid mismatch:\n got      {got}\n expected {expected}"
    print("Test 4 passed: None (the default) raises a user error outside the stored grid")


# ============================================================================
# Test 5: Spatially gridded data honours the frequency modes as well
# ============================================================================
def test_field_variants():
    field3 = arts.SpectralGeodeticField3()
    field3.grids = (arts.AscendingGrid(f_stored),
                    arts.LatGrid([-90.0, 90.0]),
                    arts.LonGrid([-180.0, 179.0]))
    # ramp along frequency, flat spatially
    field3.data = np.full((2, 2, 2), 0.0) + np.array([r_lo, r_hi]).reshape(2, 1, 1)

    expected_nearest = np.array([r_lo, ramp(5e11), r_hi])
    expected_zero = np.array([0.0, ramp(5e11), 0.0])

    for mode, expected in [(IE.Nearest, expected_nearest), (IE.Zero, expected_zero)]:
        sc = arts.LambertianSurfaceScatterer(field3)
        set_frequency_extrapolation(sc.reflectivity, mode)
        got = lambertian_brdf_i(sc, f_query)
        assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
            f"Field {mode} mismatch:\n got      {got}\n expected {expected}"

    sc = arts.LambertianSurfaceScatterer(field3)
    set_frequency_extrapolation(sc.reflectivity, IE.None_)
    try:
        lambertian_brdf_i(sc, f_query)
    except RuntimeError as e:
        assert NONE_ERROR in str(e), f"Unexpected error message: {e}"
    else:
        raise AssertionError("Field None mode must raise outside the grid")

    print("Test 5 passed: spatially gridded data honours Nearest/Zero/None")


# ============================================================================
# Test 6: Fresnel -- Nearest equals edge evaluation, None/Zero are errors
# ============================================================================
def test_fresnel_modes():
    n_lo, n_hi = 1.5, 2.5

    def specular_i(sc, f_grid):
        props = sc.get_surface_scattering_model_properties(
            surf_pt, 0.0, 0.0, f_grid, za, aa, za, aa)
        return np.array(props.brdf_matrix_specular)[:, 0, 0, 0, 0, 0, 0]

    def fresnel(mode):
        sc = arts.FresnelSurfaceScatterer(spectrum([n_lo, n_hi]))
        set_frequency_extrapolation(sc.refractive_index, mode)
        return sc

    # Nearest at out-of-grid frequencies must equal Linear evaluated at the
    # clamped (edge) frequencies -- no Fresnel formula needed, just geometry
    sc_nearest = fresnel(IE.Nearest)
    got = specular_i(sc_nearest, f_query)

    sc_linear = fresnel(IE.Linear)
    expected = specular_i(sc_linear, arts.Vector([f_stored[0], 5e11, f_stored[1]]))

    assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
        f"Fresnel Nearest != edge evaluation:\n got      {got}\n expected {expected}"

    sc_none = fresnel(IE.None_)
    try:
        specular_i(sc_none, f_query)
    except RuntimeError as e:
        assert NONE_ERROR in str(e), f"Unexpected error message: {e}"
    else:
        raise AssertionError("Fresnel None mode must raise outside the grid")

    # Zero mode gives n2 = 0 outside the grid -- a user error, the refractive
    # index must stay positive
    sc_zero = fresnel(IE.Zero)
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
    sc = arts.FlatScalarSurfaceScatterer(spectrum([r_lo, r_hi]))
    set_frequency_extrapolation(sc.reflectivity, IE.Zero)
    props = sc.get_surface_scattering_model_properties(
        surf_pt, 0.0, 0.0, f_query, za, aa, za, aa)
    got = np.array(props.brdf_matrix_specular)[:, 0, 0, 0, 0, 0, 0]
    expected = np.array([0.0, ramp(5e11), 0.0])
    assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
        f"FlatScalar Zero mismatch:\n got      {got}\n expected {expected}"
    print("Test 7 passed: FlatScalar Zero gives zero reflectivity outside the grid")


# ============================================================================
# Test 8: The lower and upper sides are independent
# ============================================================================
def test_independent_sides():
    sc = arts.LambertianSurfaceScatterer(spectrum([r_lo, r_hi]))
    sc.reflectivity.frq_low = IE.Zero
    sc.reflectivity.frq_upp = IE.Nearest
    got = lambertian_brdf_i(sc, f_query)
    expected = np.array([0.0, ramp(5e11), r_hi])
    assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
        f"Zero below / Nearest above mismatch:\n got      {got}\n expected {expected}"

    sc.reflectivity.frq_low = IE.Linear
    sc.reflectivity.frq_upp = IE.None_
    below = arts.Vector([1e9, 5e11])
    got = lambertian_brdf_i(sc, below)
    expected = np.array([ramp(f) for f in below])
    assert np.allclose(got, expected, rtol=1e-12, atol=0.0), \
        f"Linear below mismatch:\n got      {got}\n expected {expected}"
    try:
        lambertian_brdf_i(sc, f_query)
    except RuntimeError as e:
        assert NONE_ERROR in str(e), f"Unexpected error message: {e}"
    else:
        raise AssertionError("None above the grid must raise")
    print("Test 8 passed: frq_low and frq_upp act independently")


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
    test_independent_sides()
    print("\nAll extrapolation tests passed!")
