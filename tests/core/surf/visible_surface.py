"""Test visible_coordinates() for four representative terrain geometries.

All tests use a tiny sphere (a = b = 1e5 m).  Grids are centred on
(lat=0, lon=0) with spacing D = 1e-5 rad.

Visibility rule: a candidate (lat, lon) is returned by visible_coordinates()
iff
  1. It is on the same hemisphere as the observer (dot product of their
     ellipsoid outward normals > 0; only relevant for grids spanning > 90°).
  2. No intermediate terrain grid cell in the DDA ray-march lies strictly
     above the straight-line chord between observer and candidate.
  3. It is not in the 3×3 neighbourhood of the observer grid cell.

Geometry descriptions
---------------------
hill    Gaussian bump.  From the peak (highest point) every other point is
        on a downward slope.  The convex hill surface is ABOVE the chord
        from peak to any other point → DDA always blocks → 0 visible.

hole    Gaussian dip.  From the centre (lowest point) the intermediate
        terrain descends toward the observer, so the chord to moderate-
        distance rim points is always ABOVE all intermediate terrain →
        the centre can see those points.

valley  Two rectangular pits joined by a deep diagonal trench.
        Pit-A bottom → Pit-B bottom: the trench floor is lower than the
        chord → not blocked.
        Flat ground → Pit-B bottom: flat terrain (h=0) lies ABOVE the
        chord to the pit (h < 0) → blocked.

plateau Elevated flat square surrounded by flat low ground.
        Wall face (south edge, h=H_top) sees southern low ground (h=0).
        Low ground just south sees the wall.
        Plateau top centre cannot see the far low-ground rows (they are
        hidden by the step edge).
"""

import pyarts3 as pyarts
import numpy as np

# ---------------------------------------------------------------------------
# Shared grid and ellipsoid
# ---------------------------------------------------------------------------
N = 15                                    # grid points per axis
D = 1e-5                                  # grid spacing [rad]
lat = np.linspace(-(N // 2) * D, (N // 2) * D, N)
lon = np.linspace(-(N // 2) * D, (N // 2) * D, N)
ell = [1e5, 1e5]                          # tiny sphere

LAT, LON = np.meshgrid(lat, lon, indexing="ij")


def make_field(name, data):
    return pyarts.arts.GeodeticField2(
        name=name,
        grid_names=["Latitude", "Longitude"],
        grids=[lat, lon],
        data=data.tolist(),
    )


def visible_set(field, obs_lat, obs_lon):
    """Return frozenset of (lat, lon) tuples visible from the observer."""
    return frozenset(
        (round(float(c[0]), 9), round(float(c[1]), 9))
        for c in pyarts.arts.geodetic.visible_coordinates(
            pos=[obs_lat, obs_lon], ellipsoid=ell, hfield=field
        )
    )


def nearest(grid, v):
    return int(np.argmin(np.abs(grid - v)))


# ---------------------------------------------------------------------------
# 1.  HILL – Gaussian bump, sigma = 3 D, peak H = 100 m
# ---------------------------------------------------------------------------
R2 = (LAT ** 2 + LON ** 2) / D ** 2
hill_data = 100.0 * np.exp(-R2 / (2 * 3.0 ** 2))
hill_field = make_field("hill", hill_data)

peak_vis = visible_set(hill_field, 0.0, 0.0)
assert len(peak_vis) == 0, (
    f"Hill: peak should see 0 other points, got {len(peak_vis)}"
)
print(f"PASS: hill – peak sees 0 points")

# ---------------------------------------------------------------------------
# 2.  HOLE – Gaussian dip, sigma = 3 D, H = 100 m
#     h(r) = H * (exp(r^2 / 2sigma^2) - 1)  → centre = 0, rim rises outward.
# ---------------------------------------------------------------------------
hole_data = 100.0 * (np.exp(R2 / (2 * 3.0 ** 2)) - 1.0)
hole_field = make_field("hole", hole_data)

centre_vis = visible_set(hole_field, 0.0, 0.0)
assert len(centre_vis) > 0, "Hole: centre must see at least some points"

# The four cardinal points at r = 2D (h ≈ 25 m) must be visible from centre.
for la, lo in [(lat[N // 2 + 2], lat[N // 2]),
               (lat[N // 2 - 2], lat[N // 2]),
               (lat[N // 2], lat[N // 2 + 2]),
               (lat[N // 2], lat[N // 2 - 2])]:
    key = (round(la, 9), round(lo, 9))
    assert key in centre_vis, (
        f"Hole: centre must see r=2D point {key}, but does not.\n"
        f"Centre sees {len(centre_vis)} points total."
    )

print(f"PASS: hole – centre sees {len(centre_vis)} points; "
      "all four r=2D cardinal neighbours confirmed visible")

# ---------------------------------------------------------------------------
# 3.  VALLEY – two rectangular pits connected by a deep diagonal trench
#
#     Pit-A  : rows  0..3,  cols  0..3   h = -H_pit   (lower-left)
#     Pit-B  : rows 11..14, cols 11..14  h = -H_pit   (upper-right)
#     Trench : |row - col| <= 1          h = -H_trench (deeper)
#     Flat   : everywhere else           h = 0
# ---------------------------------------------------------------------------
H_pit = 300.0
H_trench = 400.0   # trench floor below pit floors

valley_data = np.zeros((N, N))
valley_data[:4, :4] = -H_pit
valley_data[11:, 11:] = -H_pit
for i in range(N):
    for j in range(N):
        if abs(i - j) <= 1:
            valley_data[i, j] = -H_trench

valley_field = make_field("valley", valley_data)

# pit-A: (row 1, col 1) — pure pit cell, not on trench diagonal
pitA_lat, pitA_lon = lat[1], lon[1]
# pit-B: (row 12, col 13) — pure pit cell
pitB_lat, pitB_lon = lat[12], lon[13]
# flat observer: (row 0, col 7)  h = 0
flat_lat, flat_lon = lat[0], lon[7]

pitA_vis = visible_set(valley_field, pitA_lat, pitA_lon)
flat_vis = visible_set(valley_field, flat_lat,  flat_lon)
pitB_key = (round(pitB_lat, 9), round(pitB_lon, 9))

assert pitB_key in pitA_vis, (
    f"Valley: pit-A ({pitA_lat:.2e},{pitA_lon:.2e}) must see "
    f"pit-B ({pitB_lat:.2e},{pitB_lon:.2e}).\n"
    f"Visible from pit-A (first 5): {sorted(pitA_vis)[:5]}"
)
assert pitB_key not in flat_vis, (
    f"Valley: flat observer ({flat_lat:.2e},{flat_lon:.2e}) must NOT see "
    f"pit-B ({pitB_lat:.2e},{pitB_lon:.2e})."
)

print(f"PASS: valley – pit-A sees pit-B; flat ground does not")

# ---------------------------------------------------------------------------
# 4.  PLATEAU – elevated 9×9 square (P_half=4) on a 15×15 flat base
#
#     Plateau : rows [i_low..i_high], cols [i_low..i_high]  h = H_top = 200 m
#     Base    : everything else                              h = 0
#     i_low = 3, i_high = 11  (grid centre ± 4 rows/cols)
# ---------------------------------------------------------------------------
P_half = 4
H_top = 200.0
i_low = N // 2 - P_half   # = 3
i_high = N // 2 + P_half   # = 11

plateau_data = np.zeros((N, N))
plateau_data[i_low:i_high + 1, i_low:i_high + 1] = H_top
plateau_field = make_field("plateau", plateau_data)

wall_lat = lat[i_low]           # south wall edge (h = H_top)
wall_lon = 0.0
low2_lat = lat[i_low - 2]       # 2 rows south of wall (h = 0, row 1)

wall_vis = visible_set(plateau_field, wall_lat, wall_lon)
low2_vis = visible_set(plateau_field, low2_lat, 0.0)
top_vis = visible_set(plateau_field,  0.0,      0.0)

# (A) wall sees some southern low ground (h = 0, lat < wall_lat)
assert any(la < wall_lat for la, _ in wall_vis), (
    f"Plateau: wall ({wall_lat:.2e},0) must see southern low ground.\n"
    f"Visible lats: {sorted(set(la for la, _ in wall_vis))}"
)

# (B) low ground (2 rows south) sees the wall face
wall_lat_r = round(wall_lat, 9)
assert any(round(la, 9) == wall_lat_r for la, _ in low2_vis), (
    f"Plateau: low-ground ({low2_lat:.2e},0) must see wall at {wall_lat:.2e}.\n"
    f"Visible lats: {sorted(set(la for la, _ in low2_vis))}"
)

# (C) plateau top cannot see the far low-ground rows
far_s = round(lat[0],  9)
far_n = round(lat[-1], 9)
top_lats = {round(la, 9) for la, _ in top_vis}
assert far_s not in top_lats and far_n not in top_lats, (
    "Plateau: top must NOT see far low-ground rows.\n"
    f"Visible lats: {sorted(top_lats)}"
)

print("PASS: plateau – wall sees low ground; low ground sees wall; "
      "top does not see far low ground")

print("\nAll visible_surface tests passed.")
