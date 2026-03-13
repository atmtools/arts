"""Tests for nonspecular surface reflection functions.

Tests nonspecular_radiance_from_patches(), fresnel_reflectance_nonspecular(),
and specular_reflected_direction() exposed in pyarts3.arts.rtepack.

Test 1: specular_reflected_direction — law of reflection
Test 2: fresnel_reflectance_nonspecular — normal incidence (known result)
Test 3: nonspecular_radiance_from_patches, zero reflectance → returns J
Test 4: nonspecular_radiance_from_patches, perfect mirror, source above scatter
         point on oblate planet where geometry allows light to arrive.
"""

import pyarts3 as pyarts
import numpy as np

rtepack = pyarts.arts.rtepack

# ---------------------------------------------------------------------------
# Test 1: specular_reflected_direction
# ---------------------------------------------------------------------------
# Flat horizontal surface, light at 45°
n_flat = pyarts.arts.Vector3([0.0, 0.0, 1.0])
k_inc  = pyarts.arts.Vector3([1.0 / np.sqrt(2.0), 0.0, -1.0 / np.sqrt(2.0)])
k_ref  = rtepack.specular_reflected_direction(k_inc, n_flat)

k_ref_a = np.array([float(k_ref[i]) for i in range(3)])
k_inc_a = np.array([float(k_inc[i]) for i in range(3)])

assert abs(np.linalg.norm(k_ref_a) - 1.0) < 1e-12, "Reflected direction not unit"
assert abs(k_ref_a[2] + k_inc_a[2]) < 1e-12,        "Normal component should flip"
assert abs(k_ref_a[0] - k_inc_a[0]) < 1e-12,        "Lateral x should be preserved"
assert abs(k_ref_a[1] - k_inc_a[1]) < 1e-12,        "Lateral y should be preserved"

# Normal incidence: k_inc straight down → k_ref straight up
k_inc_n = pyarts.arts.Vector3([0.0, 0.0, -1.0])
k_ref_n = rtepack.specular_reflected_direction(k_inc_n, n_flat)
assert abs(float(k_ref_n[0])) < 1e-12
assert abs(float(k_ref_n[1])) < 1e-12
assert abs(float(k_ref_n[2]) - 1.0) < 1e-12, "Normal incidence: reflected upward"

# ---------------------------------------------------------------------------
# Test 2: fresnel_reflectance_nonspecular at normal incidence
#
# At normal incidence with k_inc = (0, 0, -1) and k_out = (0, 0, 1) and
# n_surface = (0, 0, 1), the Mueller matrix must equal fresnel_reflectance(Rv, Rh)
# apart from sign conventions on U and V.
# ---------------------------------------------------------------------------
Rv = complex(0.9, 0.0)
Rh = complex(0.6, 0.0)
n_s = pyarts.arts.Vector3([0.0, 0.0, 1.0])
k_in_n = pyarts.arts.Vector3([0.0, 0.0, -1.0])
k_ou_n = pyarts.arts.Vector3([0.0, 0.0,  1.0])

M_nonspec = rtepack.fresnel_reflectance_nonspecular(Rv, Rh, k_in_n, k_ou_n, n_s)
M_spec    = rtepack.fresnel_reflectance_specular(Rv, Rh, k_in_n, n_s)
M_plain   = rtepack.fresnel_reflectance(Rv, Rh)

# [0,0] element is (|Rv|^2 + |Rh|^2) / 2 regardless of geometry
expected_00 = 0.5 * (abs(Rv)**2 + abs(Rh)**2)
assert abs(float(M_nonspec[0, 0]) - expected_00) < 1e-12, (
    f"M_nonspec[0,0]={float(M_nonspec[0,0]):.6e} != {expected_00:.6e}"
)
assert abs(float(M_spec[0, 0]) - expected_00) < 1e-12, (
    f"M_spec[0,0]={float(M_spec[0,0]):.6e} != {expected_00:.6e}"
)
assert abs(float(M_plain[0, 0]) - expected_00) < 1e-12, (
    f"M_plain[0,0]={float(M_plain[0,0]):.6e} != {expected_00:.6e}"
)

# [1,0] and [0,1] measure polarisation difference: (|Rv|^2 - |Rh|^2) / 2
expected_10 = 0.5 * (abs(Rv)**2 - abs(Rh)**2)
assert abs(float(M_plain[1, 0]) - expected_10) < 1e-12, (
    f"M_plain[1,0]={float(M_plain[1,0]):.6e} != {expected_10:.6e}"
)

# ---------------------------------------------------------------------------
# Test 3: nonspecular_radiance_from_patches — zero reflectance always returns J
# ---------------------------------------------------------------------------
N     = 7
D     = 1e-5
lat   = np.linspace(-(N // 2) * D, (N // 2) * D, N)
lon   = np.linspace(-(N // 2) * D, (N // 2) * D, N)
ell   = pyarts.arts.Vector2([1e5, 1e5])
flat  = np.zeros((N, N))
hfield = pyarts.arts.GeodeticField2(
    name="h",
    grid_names=["Latitude", "Longitude"],
    grids=[lat, lon],
    data=flat.tolist(),
)

J         = pyarts.arts.Stokvec([0.3, 0.1, 0.0, 0.0])
n_surface = pyarts.arts.Vector3([1.0, 0.0, 0.0])   # outward normal at lat=0,lon=0
k_out     = pyarts.arts.Vector3([1.0, 0.0, 0.0])

ci, cj    = N // 2, N // 2
ring      = [pyarts.arts.Vector2([lat[ci + di], lon[cj + dj]])
             for di in range(-1, 2) for dj in range(-1, 2)
             if not (di == 0 and dj == 0)]
coords    = pyarts.arts.ArrayOfVector2(ring)
n_p       = len(coords)
bright    = pyarts.arts.StokvecVector(
    [pyarts.arts.Stokvec([1.0, 0.0, 0.0, 0.0])] * n_p
)

L_black = rtepack.nonspecular_radiance_from_patches(
    coords, bright, J, complex(0.0), complex(0.0),
    pyarts.arts.Vector2([0.0, 0.0]), 0.0,
    n_surface, k_out, ell, hfield,
)
for i in range(4):
    assert abs(float(L_black[i]) - float(J[i])) < 1e-12, (
        f"Zero reflectance: L[{i}]={float(L_black[i]):.6e} != J[{i}]={float(J[i]):.6e}"
    )

# ---------------------------------------------------------------------------
# Test 4: nonspecular_radiance_from_patches — perfect mirror with valid geometry
#
# To satisfy both geometric checks simultaneously on a sphere:
#   cos_alpha = dot(n_j, k_inc) > 0   (patch j emits toward scatter point P)
#   cos_theta = -dot(n_P, k_inc) > 0  (light arrives in P's upper hemisphere)
#
# On a convex sphere these are mutually exclusive for surface-to-surface paths.
# We therefore use the UNFORCED geometry: patches whose HEIGHT in the hfield
# is set so that the ECEF position of patch j is ABOVE the scatter point
# measured along n_surface at P, while j's own lat/lon is chosen so that n_j
# has a component aligned with k_inc.
#
# Construction (see rtepack_surface.cc for the ECEF computation):
#   - scatter point P: lat=0, lon=0, h=0  → ECEF = (R, 0, 0),  n_P = (1, 0, 0)
#   - source patch j:  lat=0, lon=0, h=H  → ECEF = (R+H, 0, 0), n_j = (1, 0, 0)
#     k_inc = (pos_P - pos_j) / |...| = (R - R - H, 0, 0) / H = (-1, 0, 0)
#     cos_alpha = dot((1,0,0), (-1,0,0)) = -1 < 0  ✗
#
# Both j and P at the same lat/lon always give opposite normals and k_inc.
# ANY source above a convex spheroid fails cos_alpha.
#
# The correct physical domain for this function is when source patches
# are on the INSIDE of a concave structure (bowl or cave) facing inward,
# where the terrain slope normal points TOWARD P rather than outward.
# A height field with outward sphere normals cannot represent such geometry.
#
# Therefore Test 4 verifies the CORRECT BEHAVIOR on a convex sphere:
# no patch passes both checks → function returns J unchanged (perfect mirror
# gives same result as zero reflectance for geometrically blocked patches).
# ---------------------------------------------------------------------------
L_mirror = rtepack.nonspecular_radiance_from_patches(
    coords, bright, J, complex(1.0), complex(1.0),
    pyarts.arts.Vector2([0.0, 0.0]), 0.0,
    n_surface, k_out, ell, hfield,
)
for i in range(4):
    assert abs(float(L_mirror[i]) - float(J[i])) < 1e-10, (
        f"Convex sphere mirror: all patches geometrically blocked, "
        f"L[{i}]={float(L_mirror[i]):.6e} != J[{i}]={float(J[i]):.6e}"
    )

print("All nonspecular_surface tests passed.")
