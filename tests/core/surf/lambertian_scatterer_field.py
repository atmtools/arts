"""Tests for LambertianSurfaceScattererField.

Verifies:
1. Construction and attribute access.
2. At a homogeneous field (same reflectivity everywhere), the spatial variant
   produces identical results to LambertianSurfaceScatterer at any lat/lon.
3. Spatial interpolation: reflectivity varies linearly across latitude; the
   result at the midpoint matches the average of the two endpoint values.
4. XML round-trip: write and re-read the object, then verify BRDF is unchanged.
5. MapOfSurfaceScatteringModel accepts the new type.
"""

import numpy as np
import pyarts3 as pyarts
import tempfile
import os

arts = pyarts.arts

# ---------------------------------------------------------------------------
# Common grids
# ---------------------------------------------------------------------------
f_grid    = arts.Vector([1e11, 2e11, 3e11])
za_inc    = arts.Vector([0.0, 90.0])
aa_inc    = arts.Vector([0.0])
za_scat   = arts.Vector([0.0, 90.0])
aa_scat   = arts.Vector([0.0])
surf_pt   = arts.SurfacePoint()

# ---------------------------------------------------------------------------
# Test 1: Basic construction and attribute access
# ---------------------------------------------------------------------------
field3 = arts.SortedGriddedField3(
    name="reflectivity",
    grid_names=["Latitude", "Longitude", "Frequency"],
    grids=[[-90.0, 90.0], [-180.0, 180.0], [1e10, 1e12]],
    data=np.full((2, 2, 2), 0.5).tolist(),
)
scf = arts.LambertianSurfaceScattererField(field3)
assert scf.reflectivity_field.gridnames[0] == "Latitude", "Grid name mismatch"
print("Test 1 passed: construction and attribute access")


# ---------------------------------------------------------------------------
# Test 2: Homogeneous field equals uniform LambertianSurfaceScatterer
# ---------------------------------------------------------------------------
r_uniform = 0.3

spec1d = arts.SortedGriddedField1(
    name="spectrum",
    grid_names=["Frequency"],
    grids=[[1e10, 1e12]],
    data=[r_uniform, r_uniform],
)
sc1 = arts.LambertianSurfaceScatterer(spec1d)

field3_uniform = arts.SortedGriddedField3(
    name="reflectivity",
    grid_names=["Latitude", "Longitude", "Frequency"],
    grids=[[-90.0, 90.0], [-180.0, 180.0], [1e10, 1e12]],
    data=np.full((2, 2, 2), r_uniform).tolist(),
)
sc_field = arts.LambertianSurfaceScattererField(field3_uniform)

# Use lat=30, lon=45 (an interior point — should extrapolate uniformly)
props1 = sc1.get_surface_scattering_model_properties(
    surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
props_f = sc_field.get_surface_scattering_model_properties(
    surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)

brdf1  = np.array(props1.brdf_matrix)
brdf_f = np.array(props_f.brdf_matrix)
emiss1  = np.array(props1.emissivity_vector)
emiss_f = np.array(props_f.emissivity_vector)

assert np.allclose(brdf1, brdf_f, atol=1e-12), (
    f"BRDF mismatch between uniform and homogeneous field:\n"
    f"max diff = {np.max(np.abs(brdf1 - brdf_f)):.2e}")
assert np.allclose(emiss1, emiss_f, atol=1e-12), (
    f"Emissivity mismatch:\n max diff = {np.max(np.abs(emiss1 - emiss_f)):.2e}")
print("Test 2 passed: homogeneous field equals uniform scatterer")


# ---------------------------------------------------------------------------
# Test 3: Linear spatial interpolation in latitude
#
# Build a field where reflectivity varies linearly from 0.2 (lat=-90) to
# 0.8 (lat=+90), constant in longitude and frequency.  The midpoint (lat=0)
# should give exactly 0.5.
# ---------------------------------------------------------------------------
r_south = 0.2
r_north = 0.8
data_lin = np.zeros((2, 2, 2))
data_lin[0, :, :] = r_south   # lat = -90
data_lin[1, :, :] = r_north   # lat = +90

field3_lin = arts.SortedGriddedField3(
    name="reflectivity",
    grid_names=["Latitude", "Longitude", "Frequency"],
    grids=[[-90.0, 90.0], [-180.0, 180.0], [1e10, 1e12]],
    data=data_lin.tolist(),
)
sc_lin = arts.LambertianSurfaceScattererField(field3_lin)

# At lat=0, expected r = 0.5 everywhere
import math
expected_brdf = 0.5
expected_emiss = 0.5

props_mid = sc_lin.get_surface_scattering_model_properties(
    surf_pt, 0.0, 0.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
brdf_mid  = np.array(props_mid.brdf_matrix)
emiss_mid = np.array(props_mid.emissivity_vector)

# Check BRDF [f, zi, ai, zs, as, 0, 0] = r
for fi in range(len(f_grid)):
    assert abs(float(brdf_mid[fi, 0, 0, 0, 0, 0, 0]) - expected_brdf) < 1e-12, (
        f"BRDF at lat=0 f[{fi}]: {float(brdf_mid[fi,0,0,0,0,0,0]):.6e} != {expected_brdf:.6e}")
    assert abs(float(emiss_mid[fi, 0, 0]) - expected_emiss) < 1e-12, (
        f"Emissivity at lat=0 f[{fi}]: {float(emiss_mid[fi,0,0]):.6e} != {expected_emiss:.6e}")

# Also check at lat=-90 (should give r_south) and lat=+90 (r_north)
for lat_val, r_expected in [(-90.0, r_south), (90.0, r_north)]:
    props_e = sc_lin.get_surface_scattering_model_properties(
        surf_pt, lat_val, 0.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
    brdf_e  = np.array(props_e.brdf_matrix)
    for fi in range(len(f_grid)):
        val = float(brdf_e[fi, 0, 0, 0, 0, 0, 0])
        exp = r_expected
        assert abs(val - exp) < 1e-12, (
            f"BRDF at lat={lat_val}: {val:.6e} != {exp:.6e}")
print("Test 3 passed: linear spatial interpolation")


# ---------------------------------------------------------------------------
# Test 4: XML round-trip
# ---------------------------------------------------------------------------
with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
    fname = f.name
try:
    sc_field.savexml(fname)
    sc_reload = arts.LambertianSurfaceScattererField()
    sc_reload.readxml(fname)

    props_orig   = sc_field.get_surface_scattering_model_properties(
        surf_pt, 30.0, 60.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
    props_reload = sc_reload.get_bulk_surface_scattering_properties(
        surf_pt, 30.0, 60.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)

    brdf_orig   = np.array(props_orig.brdf_matrix)
    brdf_reload = np.array(props_reload.brdf_matrix)
    assert np.allclose(brdf_orig, brdf_reload, atol=1e-12), (
        f"XML round-trip: BRDF mismatch, max diff = "
        f"{np.max(np.abs(brdf_orig - brdf_reload)):.2e}")
    print("Test 4 passed: XML round-trip")
finally:
    os.unlink(fname)


# ---------------------------------------------------------------------------
# Test 5: MapOfSurfaceScatteringModel accepts LambertianSurfaceScattererField
# ---------------------------------------------------------------------------
mosm = arts.MapOfSurfaceScatteringModel()
mosm.add("albedo_field", sc_field)
assert "albedo_field" in mosm, "Model not in map"

bulk = mosm.get_surface_scattering_model_properties(
    surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
brdf_bulk = np.array(bulk.brdf_matrix)
assert not np.all(brdf_bulk == 0.0), "Bulk BRDF should be non-zero"
print("Test 5 passed: MapOfSurfaceScatteringModel accepts new type")


print("\nAll LambertianSurfaceScattererField tests passed.")
