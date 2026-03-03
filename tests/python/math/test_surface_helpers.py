import numpy as np

from pyarts3 import arts


def _expected_fresnel_matrix(Rv: complex, Rh: complex) -> np.ndarray:
    rv = np.abs(Rv) ** 2
    rh = np.abs(Rh) ** 2
    rmean = 0.5 * (rv + rh)
    rdiff = 0.5 * (rv - rh)
    a = Rh * np.conjugate(Rv)
    b = Rv * np.conjugate(Rh)
    c = 0.5 * np.real(a + b)
    d = 0.5 * np.imag(a - b)

    expected = np.zeros((4, 4), dtype=np.float64)
    expected[0, 0] = rmean
    expected[1, 0] = rdiff
    expected[0, 1] = rdiff
    expected[1, 1] = rmean
    expected[2, 2] = c
    expected[2, 3] = d
    expected[3, 2] = -d
    expected[3, 3] = c
    return expected


def _to_array(value) -> np.ndarray:
    return np.asarray(value)


results = {}

n_surface = [0.0, 0.0, 1.0]

Rv = 0.3 + 0.2j
Rh = -0.1 + 0.7j

actual = _to_array(arts.rtepack.fresnel_reflectance(Rv, Rh))
expected = _expected_fresnel_matrix(Rv, Rh)

print("Fresnel reflectance matrix (expected):\n", expected)
print("Fresnel reflectance matrix (actual):\n", actual)

assert np.allclose(actual, expected, rtol=1e-12, atol=1e-15)
results["fresnel_matrix"] = (actual, expected)


straight_in = [0.0, 0.0, -1.0]
straight_out = _to_array(
    arts.rtepack.specular_reflected_direction(straight_in, n_surface)
)
print("Specular reflection (normal incidence) computed:", straight_out)
print("Specular reflection (normal incidence) expected: [0.0, 0.0, 1.0]")
assert np.allclose(straight_out, [0.0, 0.0, 1.0])
results["specular_normal"] = straight_out

angle = np.pi / 6
oblique_in = [np.sin(angle), 0.0, -np.cos(angle)]
oblique_out = _to_array(
    arts.rtepack.specular_reflected_direction(oblique_in, n_surface)
)
expected_oblique = [np.sin(angle), 0.0, np.cos(angle)]
print("Specular reflection (oblique incidence) computed:", oblique_out)
print("Specular reflection (oblique incidence) expected:", expected_oblique)
assert np.allclose(oblique_out, expected_oblique, rtol=1e-12)
results["specular_oblique"] = oblique_out

Rv = 0.4 + 0.3j
Rh = 0.1 + 0.5j
k_inc = [0.0, 0.0, -1.0]

specular = _to_array(
    arts.rtepack.fresnel_reflectance_specular(Rv, Rh, k_inc, n_surface)
)
fresnel = _expected_fresnel_matrix(Rv, Rh)

print("Specular Fresnel rows (computed):\n", specular)
print("Specular Fresnel rows (expected with UV flip):\n", fresnel)
assert np.allclose(specular[0], fresnel[0], rtol=1e-12, atol=1e-15)
assert np.allclose(specular[1], fresnel[1], rtol=1e-12, atol=1e-15)
assert np.allclose(specular[2], -fresnel[2], rtol=1e-12, atol=1e-15)
assert np.allclose(specular[3], -fresnel[3], rtol=1e-12, atol=1e-15)
results["specular_uv_flip"] = specular

Rv = 0.2 + 0.4j
Rh = -0.3 + 0.1j
angle = np.pi / 4
k_inc = [np.sin(angle), 0.0, -np.cos(angle)]


k_out = arts.rtepack.specular_reflected_direction(k_inc, n_surface)

specular = _to_array(
    arts.rtepack.fresnel_reflectance_specular(Rv, Rh, k_inc, n_surface)
)
nonspecular = _to_array(
    arts.rtepack.fresnel_reflectance_nonspecular(
        Rv, Rh, k_inc, k_out, n_surface
    )
)

print("Specular Fresnel Mueller (computed):\n", specular)
print("Nonspecular Fresnel Mueller (computed):\n", nonspecular)
assert np.allclose(specular, nonspecular, rtol=1e-12, atol=1e-15)
results["nonspecular_agreement"] = specular

