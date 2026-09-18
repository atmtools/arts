"""Tests for MapOfSurfaceScatteringModel weighting options.

Verifies:
1. The ``Weighting`` enum is exposed to Python and its members work.
2. ``weighting_option`` defaults to ``Maximum`` and can be set / read back.
3. Behaviour: with two named models whose presence/weight differs per surface
   point, ``Maximum`` and ``Average`` produce the respective weighted blends of
   the per-model BRDF / emissivity.
"""

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts

Mosm = arts.MapOfSurfaceScatteringModel
Weighting = Mosm.Weighting

# ---------------------------------------------------------------------------
# Common grids (kept identical to the per-model evaluation domain)
# ---------------------------------------------------------------------------
f_grid = arts.Vector([1e11, 2e11, 3e11])
za_inc = arts.Vector([0.0, 90.0])
aa_inc = arts.Vector([0.0])
za_scat = arts.Vector([0.0, 90.0])
aa_scat = arts.Vector([0.0])


def lambertian(reflectivity):
    """Build a LambertianSurfaceScatterer with a constant reflectivity."""
    spectrum = arts.SortedGriddedField1(
        name="reflectivity",
        grid_names=["Frequency"],
        grids=[[1e11, 3e11]],
        data=[reflectivity, reflectivity],
    )
    return arts.LambertianSurfaceScatterer(spectrum)


def brdf_element(brdf_matrix):
    """Scalar (0,0) element of the first cell. dims [nf, nzi, nai, nzs, nas, 4, 4]"""
    return float(np.array(brdf_matrix)[0, 0, 0, 0, 0, 0, 0])


def emissivity_element(emissivity_vector):
    """Scalar (0,0) element of the first cell. dims [nf, nzs, nas, 4, 4]"""
    return float(np.array(emissivity_vector)[0, 0, 0, 0, 0])


# ---------------------------------------------------------------------------
# Test 1: Weighting enum is accessible and distinct
# ---------------------------------------------------------------------------
assert getattr(Weighting, "Maximum") is not None, "Weighting.Maximum missing"
assert getattr(Weighting, "Average") is not None, "Weighting.Average missing"
assert Weighting.Maximum != Weighting.Average, "Weighting members must differ"
print("Test 1 passed: Weighting enum is accessible and its members are distinct")


# ---------------------------------------------------------------------------
# Test 2: weighting_option defaults to Maximum
# ---------------------------------------------------------------------------
m = Mosm()
assert m.weighting_option == Weighting.Maximum, (
    f"default weighting_option is {m.weighting_option!r}, expected Maximum")
print("Test 2 passed: weighting_option defaults to Maximum")


# ---------------------------------------------------------------------------
# Test 3: weighting_option can be set to Average and read back
# ---------------------------------------------------------------------------
m.weighting_option = Weighting.Average
assert m.weighting_option == Weighting.Average, "failed to set Average"
m.weighting_option = Weighting.Maximum
assert m.weighting_option == Weighting.Maximum, "failed to reset Maximum"
print("Test 3 passed: weighting_option can be set and read back")


# ---------------------------------------------------------------------------
# Test 4: Maximum / Average produce the expected weighted blends
#
# Two models with distinct reflectivities.  The surface point carries
# per-model weights under the model keys:  w[a]=2.0, w[b]=1.0.
#
#   Model a: r = 0.2 -> BRDF_00 = 0.2, emissivity_00 = 0.8
#   Model b: r = 0.4 -> BRDF_00 = 0.4, emissivity_00 = 0.6
#
# Maximum  -> weights [2,1] -> argmax kept, rest zeroed -> [1, 0]
#            => result = 1 * prop_a = BRDF 0.2, emiss 0.8
#
# Average  -> weights [2,1]/3 -> [2/3, 1/3]
#            => BRDF  = 2/3*0.2 + 1/3*0.4 = 0.266666...
#               emiss = 2/3*0.8 + 1/3*0.6 = 0.733333...
# ---------------------------------------------------------------------------
weights = {"a": 2.0, "b": 1.0}
r_a, r_b = 0.2, 0.4

expected = {
    "Maximum": {"brdf": 1.0 * r_a, "emiss": 1.0 * (1.0 - r_a)},
    "Average": {"brdf": (2.0 * r_a + 1.0 * r_b) / 3.0,
                "emiss": (2.0 * (1.0 - r_a) + 1.0 * (1.0 - r_b)) / 3.0},
}


def bulk(opt):
    surf_point = arts.SurfacePoint()
    for key, w in weights.items():
        surf_point[arts.SurfacePropertyTag(key)] = w

    mosm = Mosm()
    mosm.add("a", lambertian(r_a))
    mosm.add("b", lambertian(r_b))
    mosm.weighting_option = getattr(Weighting, opt)
    props = mosm.get_surface_scattering_model_properties(
        surf_point, 0.0, 0.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
    return (brdf_element(props.brdf_matrix),
            emissivity_element(props.emissivity_vector))


for opt in ("Maximum", "Average"):
    brdf, emiss = bulk(opt)
    expected_outcomes = expected[opt]
    assert abs(brdf - expected_outcomes["brdf"]) < 1e-12, (
        f"{opt}: BRDF {brdf:.6e} != expected {expected_outcomes['brdf']:.6e}")
    assert abs(emiss - expected_outcomes["emiss"]) < 1e-12, (
        f"{opt}: emissivity {emiss:.6e} != expected {expected_outcomes['emiss']:.6e}")
    print(f"Test 4 ({opt}) passed: weighted blend matches expected value")

# The two modes must be distinguishable for this configuration.
assert bulk("Maximum") != bulk("Average"), "Maximum and Average coincided"

print("\nAll MapOfSurfaceScatteringModel weighting tests passed.")
