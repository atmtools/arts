"""Tests for FresnelSurfaceScatterer and FresnelSurfaceScattererField.

Verifies:
1. Construction and attribute access (n1 default 1.0).
2. Known values: n1=1, n2=1.333, za_inc=[0, 30] — compare the BRDF specular
   muelmat against the closed-form Fresnel coefficients recomputed here, and
   the complementary (I4 - M) emissivity.  Diffuse tensors stay zero.
3. Grazing incidence (za=90): Rv = Rh = -1 => BRDF = 1, emissivity = 0.
4. Field variant homogeneous == spectrum variant at an interior lat/lon.
5. Total internal reflection: n1=2, n2=1, za_inc=60 => Rv = Rh = 1.
6. XML round-trip.
7. MapOfSurfaceScatteringModel accepts the new types.
"""

import cmath
import math
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

# Stored spectrum (frequency) and simulation f_grid
n2_water = 1.333
spec = arts.SortedGriddedField1(
    name="n2",
    grid_names=["Frequency"],
    grids=[[1e10, 1e12]],
    data=[n2_water, n2_water],
)


def fresnel_amplitudes(n1, n2, za_deg):
    """Closed-form complex Fresnel amplitude coefficients (theta in degrees).

    Mirrors physics_funcs.cc pair overload, including the TIR guard.
    """
    th = math.radians(za_deg)
    sin2 = n1 * math.sin(th) / n2
    if abs(sin2) > 1.0:
        return 1.0, 1.0  # total internal reflection
    cos1 = math.cos(th)
    cos2 = math.cos(math.asin(sin2))
    Rv = (n2 * cos1 - n1 * cos2) / (n2 * cos1 + n1 * cos2)
    Rh = (n1 * cos1 - n2 * cos2) / (n1 * cos1 + n2 * cos2)
    return Rv, Rh


def mun_fresnel_reflectance(Rv, Rh):
    """Mirror of rtepack::fresnel_reflectance (rtepack_surface.cc)."""
    rv = abs(Rv) ** 2
    rh = abs(Rh) ** 2
    rmean = 0.5 * (rv + rh)
    rdiff = 0.5 * (rv - rh)
    c = (Rh * Rv.conjugate()).real
    d = 0.5 * (Rh * Rv.conjugate() - Rv * Rh.conjugate()).imag
    M = np.zeros((4, 4))
    M[0, 0] = rmean
    M[1, 1] = rmean
    M[0, 1] = rdiff
    M[1, 0] = rdiff
    M[2, 2] = c
    M[3, 3] = c
    M[2, 3] = d
    M[3, 2] = -d
    return M


# ---------------------------------------------------------------------------
# Test 1: Construction and attribute access
# ---------------------------------------------------------------------------
sc = arts.FresnelSurfaceScatterer(spec, 1.0)
assert sc.n1 == 1.0, "n1 mismatch"
assert np.array_equal(np.array(sc.refractive_index_spectrum.data), [n2_water, n2_water])

sc_def = arts.FresnelSurfaceScatterer(spec)
assert sc_def.n1 == 1.0, "n1 default must be 1.0"
print("Test 1 passed: construction and attribute access")

# ---------------------------------------------------------------------------
# Test 2: Known values compared to the closed form
# ---------------------------------------------------------------------------
props = sc.get_surface_scattering_model_properties(
    surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)

brdf = np.array(props.brdf_matrix_specular)
emiss = np.array(props.emissivity_vector_specular)
brdf_d = np.array(props.brdf_matrix_diffuse)
emiss_d = np.array(props.emissivity_vector_diffuse)

assert brdf.shape == (len(f_grid), len(za_inc), len(aa_inc), len(za_scat), len(aa_scat), 4, 4)
assert emiss.shape == (len(f_grid), len(za_scat), len(aa_scat), 4, 4)

assert np.all(brdf_d == 0.0), "Diffuse BRDF must be zero for the Fresnel model"
assert np.all(emiss_d == 0.0), "Diffuse emissivity must be zero for the Fresnel model"

n1 = 1.0
for fi in range(len(f_grid)):
    for zi, za in enumerate(za_inc):
        Rv, Rh = fresnel_amplitudes(n1, n2_water, za)
        M = mun_fresnel_reflectance(Rv, Rh)
        for ai in range(len(aa_inc)):
            for zs in range(len(za_scat)):
                for as_ in range(len(aa_scat)):
                    assert np.allclose(brdf[fi, zi, ai, zs, as_], M, atol=1e-12), (
                        f"BRDF mismatch at f[{fi}] za={za}:\n"
                        f"{brdf[fi, zi, ai, zs, as_]}\n!=\n{M}")
    for zs, za in enumerate(za_scat):
        Rv, Rh = fresnel_amplitudes(n1, n2_water, za)
        M = mun_fresnel_reflectance(Rv, Rh)
        E = np.eye(4) - M
        for as_ in range(len(aa_scat)):
            assert np.allclose(emiss[fi, zs, as_], E, atol=1e-12), (
                f"Emissivity mismatch at f[{fi}] za={za}:\n"
                f"{emiss[fi, zs, as_]}\n!=\n{E}")

# Emissivity is the exact complement: M + E = I4
for fi in range(len(f_grid)):
    M0 = mun_fresnel_reflectance(*fresnel_amplitudes(n1, n2_water, za_scat[0]))
    assert np.allclose(brdf[fi, 0, 0, 0, 0] + emiss[fi, 0, 0], np.eye(4), atol=1e-12), (
        "Fresnel BRDF + emissivity must equal the identity (energy consistency)")
print("Test 2 passed: known values match the closed form")

# ---------------------------------------------------------------------------
# Test 3: Grazing incidence — Rv = Rh = -1 => total reflection, emittance 0
# ---------------------------------------------------------------------------
za_graze = arts.Vector([90.0])
scg = arts.FresnelSurfaceScatterer(spec)
props_g = scg.get_surface_scattering_model_properties(
    surf_pt, 0.0, 0.0, arts.Vector([2e11]), za_graze, aa_inc,
    za_graze, aa_scat)
brdf_g = np.array(props_g.brdf_matrix_specular)
emiss_g = np.array(props_g.emissivity_vector_specular)
# Incident side: perfect reflection at grazing
assert abs(float(brdf_g[0, 0, 0, 0, 0, 0, 0]) - 1.0) < 1e-12
# Emissivity side at za_scat=90 (mirror of grazing): 1 - 1 = 0
assert abs(float(emiss_g[0, 0, 0, 0, 0])) < 1e-12, (
    "Emissivity must be 0 at grazing incidence (total reflection)")
print("Test 3 passed: grazing incidence gives total reflection")

# ---------------------------------------------------------------------------
# Test 4: Field variant homogeneous == spectrum variant
# ---------------------------------------------------------------------------
field3 = arts.SortedGriddedField3(
    name="n2",
    grid_names=["Latitude", "Longitude", "Frequency"],
    grids=[[-90.0, 90.0], [-180.0, 180.0], [1e10, 1e12]],
    data=np.full((2, 2, 2), n2_water).tolist(),
)
scf = arts.FresnelSurfaceScattererField(field3, 1.0)
assert scf.n1 == 1.0

props1 = sc.get_surface_scattering_model_properties(
    surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
props_f = scf.get_surface_scattering_model_properties(
    surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)

assert np.allclose(np.array(props1.brdf_matrix_specular),
                   np.array(props_f.brdf_matrix_specular), atol=1e-12), (
    "Field variant must equal the spectrum variant for a homogeneous field")
assert np.allclose(np.array(props1.emissivity_vector_specular),
                   np.array(props_f.emissivity_vector_specular), atol=1e-12)
print("Test 4 passed: homogeneous field equals spectrum variant")

# ---------------------------------------------------------------------------
# Test 5: Total internal reflection (n1 = 2, n2 = 1, za_inc = 60)
# ---------------------------------------------------------------------------
spec_tir = arts.SortedGriddedField1(
    name="n2",
    grid_names=["Frequency"],
    grids=[[1e10, 1e12]],
    data=[1.0, 1.0],
)
sc_tir = arts.FresnelSurfaceScatterer(spec_tir, 2.0)
assert sc_tir.n1 == 2.0
za60 = arts.Vector([60.0])
props_t = sc_tir.get_surface_scattering_model_properties(
    surf_pt, 0.0, 0.0, arts.Vector([2e11]), za60, aa_inc, za60, aa_scat)
brdf_t = np.array(props_t.brdf_matrix_specular)
emiss_t = np.array(props_t.emissivity_vector_specular)
# TIR => Rv = Rh = 1 => BRDF = identity muelmat (all-ones on the diagonal),
# emissivity = 0
assert float(brdf_t[0, 0, 0, 0, 0, 0, 0]) == 1.0
assert float(brdf_t[0, 0, 0, 0, 0, 1, 0]) == 0.0
assert abs(float(emiss_t[0, 0, 0, 0, 0])) < 1e-12
print("Test 5 passed: total internal reflection")

# ---------------------------------------------------------------------------
# Test 6: XML round-trip
# ---------------------------------------------------------------------------
with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
    fname = f.name
try:
    scf.savexml(fname)
    sc_reload = arts.FresnelSurfaceScattererField()
    sc_reload.readxml(fname)
    assert sc_reload.n1 == scf.n1
    props_r = sc_reload.get_surface_scattering_model_properties(
        surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
    assert np.allclose(np.array(props_f.brdf_matrix_specular),
                       np.array(props_r.brdf_matrix_specular), atol=1e-12), (
        "XML round-trip: BRDF mismatch")
    print("Test 6 passed: XML round-trip")
finally:
    os.unlink(fname)

# ---------------------------------------------------------------------------
# Test 7: MapOfSurfaceScatteringModel accepts the new types
# ---------------------------------------------------------------------------
mosm = arts.MapOfSurfaceScatteringModel()
mosm.add("fresnel", sc)
mosm.add("fresnel_field", scf)
assert "fresnel" in mosm and "fresnel_field" in mosm
bulk = mosm.get_surface_scattering_model_properties(
    surf_pt, 30.0, 45.0, f_grid, za_inc, aa_inc, za_scat, aa_scat)
brdf_bulk = np.array(bulk.brdf_matrix_specular)
assert not np.all(brdf_bulk == 0.0), "Bulk specular BRDF should be non-zero"
print("Test 7 passed: MapOfSurfaceScatteringModel accepts new types")

# ---------------------------------------------------------------------------
# Test 8: za > 90 must fold onto 180 - za (regression test for the raw-angle
# bug: R(180-za) = 1/R(za) if passed unfolded — e.g. za=127 blew up to ~2e6
# and za=180 gave R ~ 50 instead of ~ 0.02).
# ---------------------------------------------------------------------------
za_down = arts.Vector([91.0, 127.0, 180.0])
za_fold = arts.Vector([89.0, 53.0, 0.0])
props_d8 = sc.get_surface_scattering_model_properties(
    surf_pt, 0.0, 0.0, f_grid, za_down, aa_inc, za_down, aa_scat)
props_f8 = sc.get_surface_scattering_model_properties(
    surf_pt, 0.0, 0.0, f_grid, za_fold, aa_inc, za_fold, aa_scat)

brdf_d8 = np.array(props_d8.brdf_matrix_specular)
brdf_f8 = np.array(props_f8.brdf_matrix_specular)
emiss_d8 = np.array(props_d8.emissivity_vector_specular)
emiss_f8 = np.array(props_f8.emissivity_vector_specular)

assert np.allclose(brdf_d8, brdf_f8, atol=1e-12), (
    "BRDF at za>90 must equal BRDF at 180-za")
assert np.allclose(emiss_d8, emiss_f8, atol=1e-12), (
    "Emissivity at za>90 must equal emissivity at 180-za")

# Closed-form sanity: za=180 -> same as za=0 -> R = ((n1-n2)/(n1+n2))^2
R_normal = ((n1 - n2_water) / (n1 + n2_water)) ** 2
assert abs(float(brdf_d8[0, 2, 0, 2, 0, 0, 0]) - R_normal) < 1e-12, (
    f"za=180 BRDF[0,0] = {float(brdf_d8[0, 2, 0, 2, 0, 0, 0])}, expected {R_normal}")
assert abs(float(emiss_d8[0, 2, 0, 0, 0]) - (1.0 - R_normal)) < 1e-12, (
    "za=180 emissivity must be 1 - R_normal, never negative")
assert np.all(emiss_d8[..., 0, 0] >= 0.0), "Fresnel emissivity must never be negative"
print("Test 8 passed: za>90 folds onto 180-za")

print("\nAll FresnelSurfaceScatterer tests passed.")
