import numpy as np

from pyarts3 import arts


def _to_array(value) -> np.ndarray:
    return np.asarray(value)


results = {}

n_surface = [0.0, 0.0, 1.0]

Rv = 0.3 + 0.2j
Rh = -0.1 + 0.7j

actual = _to_array(arts.rtepack.fresnel_reflectance(Rv, Rh))
# |Rv|^2 = 0.13, |Rh|^2 = 0.5 and Rh conj(Rv) = 0.11 + 0.23j
expected = np.array(
    [
        [0.315, -0.185, 0.0, 0.0],
        [-0.185, 0.315, 0.0, 0.0],
        [0.0, 0.0, 0.11, 0.23],
        [0.0, 0.0, -0.23, 0.11],
    ]
)

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
# |Rv|^2 = 0.25, |Rh|^2 = 0.26 and Rh conj(Rv) = 0.19 + 0.17j
fresnel = np.array(
    [
        [0.255, -0.005, 0.0, 0.0],
        [-0.005, 0.255, 0.0, 0.0],
        [0.0, 0.0, 0.19, 0.17],
        [0.0, 0.0, -0.17, 0.19],
    ]
)

# On a flat surface the planes of incidence and the meridional planes coincide, so the
# reflection between the meridional bases of the two rays is the Fresnel matrix itself
print("Specular Fresnel Mueller (computed):\n", specular)
print("Specular Fresnel Mueller (expected):\n", fresnel)
assert np.allclose(specular, fresnel, rtol=1e-12, atol=1e-15)
results["specular_flat"] = specular

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

