"""Tests for FlatScalarSurfaceScatterer and FlatScalarSurfaceScattererField.

Verifies:
1. Construction and attribute access.
2. Properties: brdf_matrix_specular = Muelmat{R} at [0,0] (all angles),
   all other elements 0; emissivity_vector_specular[0,0] = 1 - R; diffuse
   tensors all zero.
3. Clamping: out-of-range stored R values are clamped to [0, 1].
4. Field variant homogeneous == spectrum variant at an interior lat/lon.
5. XML round-trip.
6. MapOfSurfaceScatteringModel accepts the new types.
"""

import os
import tempfile

import numpy as np
import pyarts3 as pyarts

arts = pyarts.arts

# ---------------------------------------------------------------------------
# Common grids
# ---------------------------------------------------------------------------
f_grid  = arts.Vector([1e11, 2e11, 3e11])
za_inc  = arts.Vector([0.0, 30.0])
aa_inc  = arts.Vector([0.0])
za_scat = arts.Vector([0.0, 90.0])
aa_scat = arts.Vector([0.0])
surf_pt = arts.SurfacePoint()

r_uniform = 0.3
spec = arts.SortedGriddedField1(
    name="reflectivity",
    grid_names=["Frequency"],
    grids=[[1e10, 1e12]],
    data=[r_uniform, r_uniform],
)


# ---------------------------------------------------------------------------
# Test 1: Construction and attribute access
# ---------------------------------------------------------------------------
sc = arts.FlatScalarSurfaceScatterer(spec)
assert np.array_equal(np.array(sc.reflectivity_spectrum.data), [r_uniform, r_uniform])
print("Test 1 passed: construction and attribute access")

# ---------------------------------------------------------------------------
# Test 2: Spectral BRDF / emissivity values
# ---------------------------------------------------------------------------
props = sc.get_surface_scattering_model_properties(
    surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)

brdf = np.array(props.brdf_matrix_specular)
emiss = np.array(props.emissivity_vector_specular)
brdf_d = np.array(props.brdf_matrix_diffuse)
emiss_d = np.array(props.emissivity_vector_diffuse)

assert brdf.shape == (len(f_grid), len(za_inc), len(aa_inc), len(za_scat), len(aa_scat), 4, 4)

# Diffuse tensors must stay zero
assert np.all(brdf_d == 0.0), "Diffuse BRDF must be zero for the flat-scalar model"
assert np.all(emiss_d == 0.0), "Diffuse emissivity must be zero for the flat-scalar model"

r_diag = np.diag([r_uniform] * 4)
e_diag = np.diag([1.0 - r_uniform] * 4)

for fi in range(len(f_grid)):
    # Specular BRDF: Muelmat{R} — full diagonal R, zero off-diagonal, all angles
    for zi in range(len(za_inc)):
        for ai in range(len(aa_inc)):
            for zs in range(len(za_scat)):
                for as_ in range(len(aa_scat)):
                    m = brdf[fi, zi, ai, zs, as_]
                    assert np.allclose(m, r_diag), f"Specular BRDF at slot must be Muelmat{{{r_uniform}}}\n{m}\n!=\n{r_diag}"
    # Specular emissivity: Muelmat{1 - R}
    for zs in range(len(za_scat)):
        for as_ in range(len(aa_scat)):
            me = emiss[fi, zs, as_]
            assert np.allclose(me, e_diag), f"Specular emissivity must be Muelmat{{{1.0 - r_uniform}}}\n{me}\n!=\n{e_diag}"

# Energy consistency: M + E = I4
for fi in range(len(f_grid)):
    assert np.allclose(
        brdf[fi, 0, 0, 0, 0] + emiss[fi, 0, 0], np.eye(4), atol=1e-15)
print("Test 2 passed: flat-scalar BRDF/emissivity values")

# ---------------------------------------------------------------------------
# Test 3: Clamping of out-of-range reflectivity
# ---------------------------------------------------------------------------
spec_clamp = arts.SortedGriddedField1(
    name="reflectivity",
    grid_names=["Frequency"],
    grids=[[1e10, 1e12]],
    data=[-1.0, 2.0],
)
sc_c = arts.FlatScalarSurfaceScatterer(spec_clamp)
props_c = sc_c.get_surface_scattering_model_properties(
    surf_pt, 0.0, 0.0, arts.Vector([1e10, 1e12]), za_inc, aa_inc, za_scat, aa_scat)
brdf_c = np.array(props_c.brdf_matrix_specular)
emiss_c = np.array(props_c.emissivity_vector_specular)
assert float(brdf_c[0, 0, 0, 0, 0, 0, 0]) == 0.0, "R=-1 must clamp to 0"
assert float(brdf_c[1, 0, 0, 0, 0, 0, 0]) == 1.0, "R=2 must clamp to 1"
assert float(emiss_c[0, 0, 0, 0, 0]) == 1.0
assert float(emiss_c[1, 0, 0, 0, 0]) == 0.0
print("Test 3 passed: clamping of out-of-range values")

# ---------------------------------------------------------------------------
# Test 4: Field variant homogeneous == spectrum variant
# ---------------------------------------------------------------------------
field3 = arts.SortedGriddedField3(
    name="reflectivity",
    grid_names=["Latitude", "Longitude", "Frequency"],
    grids=[[-90.0, 90.0], [-180.0, 180.0], [1e10, 1e12]],
    data=np.full((2, 2, 2), r_uniform).tolist(),
)
scf = arts.FlatScalarSurfaceScattererField(field3)

props1 = sc.get_surface_scattering_model_properties(
    surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
props_f = scf.get_surface_scattering_model_properties(
    surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
assert np.allclose(np.array(props1.brdf_matrix_specular),
                   np.array(props_f.brdf_matrix_specular), atol=1e-12)
assert np.allclose(np.array(props1.emissivity_vector_specular),
                   np.array(props_f.emissivity_vector_specular), atol=1e-12)
print("Test 4 passed: homogeneous field equals spectrum variant")

# ---------------------------------------------------------------------------
# Test 5: XML round-trip
# ---------------------------------------------------------------------------
with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
    fname = f.name
try:
    scf.savexml(fname)
    sc_reload = arts.FlatScalarSurfaceScattererField()
    sc_reload.readxml(fname)
    props_r = sc_reload.get_surface_scattering_model_properties(
        surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
    assert np.allclose(np.array(props_f.brdf_matrix_specular),
                       np.array(props_r.brdf_matrix_specular), atol=1e-12)
    print("Test 5 passed: XML round-trip")
finally:
    os.unlink(fname)

# ---------------------------------------------------------------------------
# Test 6: MapOfSurfaceScatteringModel accepts the new types
# ---------------------------------------------------------------------------
mosm = arts.MapOfSurfaceScatteringModel()
mosm.add("flat_scalar", sc)
mosm.add("flat_scalar_field", scf)
assert "flat_scalar" in mosm and "flat_scalar_field" in mosm
bulk = mosm.get_surface_scattering_model_properties(
    surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
brdf_bulk = np.array(bulk.brdf_matrix_specular)
assert not np.all(brdf_bulk == 0.0), "Bulk specular BRDF should be non-zero"
print("Test 6 passed: MapOfSurfaceScatteringModel accepts new types")

print("\nAll FlatScalarSurfaceScatterer tests passed.")
