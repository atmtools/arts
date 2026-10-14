#!/usr/bin/env python3

"""VDISORT with a polarized beam and a polarized surface, against field-level closed forms.

The polarized comparisons with RT3 and RT4 never excite the Q, U and V
columns of the beam, because both have an unpolarized beam, and RT3 allows
only a Lambertian surface with a beam.  This test covers both through
pyarts3 with references built from electric fields, with rtepack's
polarization primitives (meridional_basis and the Mueller matrix of a
Jones matrix, tested against fields in test_rtepack_polarization) but not VDISORT's phase
matrices or rtepack's reflection matrices:

A. Single scattering of a polarized beam in an optically thin Rayleigh
   layer.  The Rayleigh Mueller matrix comes from the dipole field
   E_out = E_in - k_out (k_out . E_in) projected on the meridional bases.
   VDISORT's diffuse operator is zero, so its formal solution is exactly the
   single-scattering integral.
B. VDISORT's raw polarized Cox-Munk reflection against the reflection of
   the fields by the specular facet: Fresnel amplitudes in the facet's
   s and p directions, projected on the meridional bases.  The scalar M00
   (the facet slope distribution and shadowing) is DISORT's, validated
   against the Fortran DISORT problems; B tests the polarization structure.
C. The beam, unpolarized and polarized, reflected by a Cox-Munk surface
   under a transparent atmosphere.  The upward radiance at the top is
   mu0 F rho S_b, Fourier-truncated as the solver's surface modes are.

Bases and Stokes vector: z up; a ray propagating along k has the
meridional basis h = k x z / |k x z| and v = h x k of
arts.rtepack.meridional_basis; I = |E_v|^2 + |E_h|^2,
Q = |E_v|^2 - |E_h|^2, U = 2 Re(E_v E_h*) and V = 2 Im(E_v E_h*).  U is pinned by the RT3
comparison.  V = 2 Im(E_v E_h*) is the convention of rtepack's Fresnel
matrices; whether it makes V positive for left-hand circular polarization is
not tested here.
"""

import numpy as np

from pyarts3 import arts

vd = arts.vdisort
rtepack = arts.rtepack


def direction(mu, phi):
    """The propagation direction of a ray at zenith cosine mu (> 0 up) and azimuth phi."""
    s = np.sqrt(max(0.0, 1.0 - mu * mu))
    return np.array([s * np.cos(phi), s * np.sin(phi), mu])


def basis(mu, phi):
    """The meridional basis (v, h) of the ray direction(mu, phi), at the poles that of the azimuth phi."""
    v, h = rtepack.meridional_basis(direction(mu, phi), np.array([np.sin(phi), -np.cos(phi), 0.0]))
    return np.array(v), np.array(h)


def mueller(jones):
    """The Mueller matrix of a 2 x 2 Jones matrix in (v, h) bases."""
    return np.array(arts.Muelmat(arts.Jones(np.asarray(jones, dtype=complex))))


def jones(out_basis, in_basis, field):
    """The Jones matrix of the linear map field(E_in) -> E_out between two bases."""
    return np.array([[a @ field(b) for b in in_basis] for a in out_basis])


def rayleigh(mu_out, phi_out, mu_in, phi_in):
    """The normalized Rayleigh phase matrix between meridional bases (M00 = 3/4 (1 + cos^2))."""
    k_in, k_out = direction(mu_in, phi_in), direction(mu_out, phi_out)
    dipole = lambda e: e - k_out * (k_out @ e)
    return 1.5 * mueller(jones(basis(mu_out, phi_out), basis(mu_in, phi_in), dipole))


def facet_reflection(mu_out, mu_in, relative_azimuth, n):
    """The Mueller reflection of the facet that reflects -mu_in (azimuth 0) into mu_out (azimuth relative_azimuth).

    The facet normal bisects the two directions.  With s perpendicular to the
    plane of incidence and p = s x k for both rays, the reflected field is
    r_s (E . s) s + r_p (E . p_in) p_out, where r_p -> +1 and r_s -> -1 for a
    perfect conductor (tangential E vanishes at the surface).
    """
    k_in, k_out = direction(-mu_in, 0.0), direction(mu_out, relative_azimuth)
    normal = (k_out - k_in) / np.linalg.norm(k_out - k_in)
    cos_local = -k_in @ normal
    s = np.cross(k_in, normal)
    s /= np.linalg.norm(s)
    p_in, p_out = np.cross(s, k_in), np.cross(s, k_out)
    root = np.sqrt(n * n - (1.0 - cos_local**2) + 0j)
    r_p = (n * n * cos_local - root) / (n * n * cos_local + root)
    r_s = (cos_local - root) / (cos_local + root)
    field = lambda e: r_s * (e @ s) * s + r_p * (e @ p_in) * p_out
    return mueller(jones(basis(mu_out, relative_azimuth), basis(-mu_in, 0.0), field))


def max_dev(a, b, scale):
    return np.abs(np.asarray(a) - np.asarray(b)).max() / scale


def test_jones_rayleigh_matches_closed_form():
    """The field construction reproduces the scattering-plane Rayleigh matrix along the principal plane."""
    for mu_out in (-0.7, 0.2, 0.9):
        Z = rayleigh(mu_out, 0.0, -0.6, 0.0)
        cos_t = direction(mu_out, 0.0) @ direction(-0.6, 0.0)
        expected = np.zeros((4, 4))
        expected[0, 0] = expected[1, 1] = 0.75 * (1 + cos_t**2)
        expected[0, 1] = expected[1, 0] = -0.75 * (1 - cos_t**2)
        expected[2, 2] = expected[3, 3] = 1.5 * cos_t
        # In the principal plane the meridional and scattering-plane bases agree up to a sign of
        # both basis vectors of a ray, which leaves the Stokes vector unchanged
        np.testing.assert_allclose(Z, expected, atol=1e-14)


def test_polarized_beam_single_scattering():
    """A. Single scattering of a polarized beam, every mode and both combined systems."""
    nquad, nfourier = 16, 3  # Rayleigh scattering has modes m <= 2 only
    depth, mu0, phi0 = 1e-6, 0.6, 0.7
    beam = 1.3 * np.array([1.0, 0.3, -0.4, 0.5])

    model0 = arts.cppvdisort(
        tau_arr=np.array([depth]),
        omega_arr=np.array([1.0]),
        NQuad=nquad,
        phase_matrix=np.zeros((2, nfourier, 1, nquad, nquad, 4, 4)),
        mu0=mu0,
        beam_stokes=np.zeros(4),
        phi0=phi0,
    )
    mu = np.asarray(model0.mu)

    # C^m and S^m of Z(mu_i, 0 <- -mu0, phi') by the trapezoid rule, exact for these trigonometric polynomials
    nphi = 64
    phis = 2 * np.pi * np.arange(nphi) / nphi
    cosine = np.zeros((nfourier, 1, nquad, 4, 4))
    sine = np.zeros_like(cosine)
    for i, m_out in enumerate(mu):
        for phi in phis:
            Z = rayleigh(m_out, 0.0, -mu0, phi)
            for m in range(nfourier):
                cosine[m, 0, i] += np.cos(m * phi) / nphi * Z
                sine[m, 0, i] += np.sin(m * phi) / nphi * Z

    def solve(stokes):
        return arts.cppvdisort(
            tau_arr=np.array([depth]),
            omega_arr=np.array([1.0]),
            NQuad=nquad,
            phase_matrix=np.zeros((2, nfourier, 1, nquad, nquad, 4, 4)),
            mu0=mu0,
            beam_stokes=stokes,
            phi0=phi0,
            beam_phase_matrix=vd.combine_beam_phase_matrices(cosine, sine),
        )

    model = solve(beam)
    phi_out = np.array([0.0, 0.3, phi0, 1.9, np.pi, 4.4])
    u = np.asarray(model.u(tau=np.array([0.0, depth]), phi=phi_out))
    n = nquad // 2

    def expected(phi_list, stokes):
        out = np.zeros((2, len(phi_list), nquad, 4))
        for k, phi in enumerate(phi_list):
            for i in range(n):
                m = mu[i]
                # Upward at the top and downward at the bottom: the attenuated beam source integrated along the ray
                up = mu0 / (mu0 + m) * -np.expm1(-depth * (1 / mu0 + 1 / m))
                down = (np.exp(-depth / mu0) - np.exp(-depth / m)) / (1 - m / mu0)
                out[0, k, i] = up / (4 * np.pi) * rayleigh(m, phi, -mu0, phi0) @ stokes
                out[1, k, n + i] = down / (4 * np.pi) * rayleigh(-m, phi, -mu0, phi0) @ stokes
        return out

    ref = expected(phi_out, beam)
    scale = np.abs(ref[..., 0]).max()
    selected = np.zeros_like(u, dtype=bool)
    selected[0, :, :n] = True  # upward at the top
    selected[1, :, n:] = True  # downward at the bottom
    deviation = np.abs(np.where(selected, u - ref, 0.0)).max(axis=(0, 1, 2)) / scale
    print(f"A. polarized beam, single scattering: max |VDISORT - closed form| / max I = {deviation}")
    assert deviation.max() < 1e-8, deviation

    # Every component is exercised, and the comparison is not blind to the beam's polarization or the azimuth sense
    polarization = np.abs(ref).max(axis=(0, 1, 2)) / scale
    assert polarization.min() > 0.05, polarization
    unpolarized = np.abs(np.where(selected, np.asarray(solve(np.array([beam[0], 0, 0, 0])).u(np.array([0.0, depth]), phi_out)) - ref, 0)).max() / scale
    mirrored = np.abs(np.where(selected, u - expected(2 * phi0 - phi_out, beam), 0)).max() / scale
    print(f"   an unpolarized beam misses by {unpolarized:.3f}, the mirrored azimuth by {mirrored:.3f}")
    assert unpolarized > 1e4 * 1e-8 and mirrored > 1e4 * 1e-8


def test_cox_munk_polarization_structure():
    """B. VDISORT's raw polarized Cox-Munk matrix against the field reflection by the facet."""
    rng = np.random.default_rng(20261013)
    for n in (1.34, 1.5 + 0.3j, 7.0 + 3.0j):
        worst = 0.0
        for _ in range(200):
            mu_out, mu_in = rng.uniform(0.1, 1.0, 2)
            azimuth = rng.uniform(0.0, 2 * np.pi)
            raw = np.asarray(vd.cox_munk_reflection(mu_out, mu_in, azimuth, 7.0, complex(n), True))
            if raw[0, 0] <= 1e-12:  # shadowed or out of the facet distribution's range
                continue
            geo = facet_reflection(mu_out, mu_in, azimuth, n)
            worst = max(worst, np.abs(raw / raw[0, 0] - geo / geo[0, 0]).max())
        print(f"B. Cox-Munk, n = {n}: max |raw / raw_00 - facet / facet_00| = {worst:.1e}")
        assert worst < 1e-12, (n, worst)


def test_beam_reflection_by_cox_munk():
    """C. The beam, unpolarized and polarized, reflected by Cox-Munk under a transparent atmosphere."""
    nquad, nfourier = 16, 24
    wind, n_index = 7.0, complex(1.5, 0.3)
    mu0, phi0, depth = 0.6, 0.4, 1e-12
    streams = np.asarray(arts.math.pdleggauss(nquad)[0])
    modes = vd.cox_munk_fourier_modes(
        wind, n_index, True, nfourier, streams, streams, mu0, azimuth_quadrature_points=2000
    )
    phi_out = np.array([0.0, 0.4, 1.2, 2.5, np.pi, 4.0, 5.5])
    n = nquad // 2

    for label, stokes in (
        ("unpolarized", 2.0 * np.array([1.0, 0.0, 0.0, 0.0])),
        ("polarized", 2.0 * np.array([1.0, -0.2, 0.5, -0.6])),
    ):
        model = arts.cppvdisort(
            tau_arr=np.array([depth]),
            omega_arr=np.array([0.0]),
            NQuad=nquad,
            phase_matrix=np.zeros((2, nfourier, 1, nquad, nquad, 4, 4)),
            mu0=mu0,
            beam_stokes=stokes,
            phi0=phi0,
            surface=modes,
        )
        mu = np.asarray(model.mu)[:n]
        u = np.asarray(model.u(tau=np.array([0.0]), phi=phi_out))[0, :, :n]

        # mu0 F rho(mu, mu0, phi - phi0) S_b, Fourier-truncated to the solver's modes in the relative azimuth
        ngrid = 1024
        grid = 2 * np.pi * np.arange(ngrid) / ngrid
        ref = np.zeros((len(phi_out), n, 4))
        mirrored = np.zeros_like(ref)  # at the mirrored relative azimuth, which must miss
        for i, m in enumerate(mu):
            g = np.array([np.asarray(vd.cox_munk_reflection(m, mu0, a, wind, n_index, True)) @ stokes for a in grid])
            coef = np.fft.rfft(g, axis=0) / ngrid
            series = lambda rel: coef[0].real + sum(
                2 * (coef[j].real * np.cos(j * rel) - coef[j].imag * np.sin(j * rel)) for j in range(1, nfourier)
            )
            for k, phi in enumerate(phi_out):
                ref[k, i] = mu0 * series(phi - phi0)
                mirrored[k, i] = mu0 * series(phi0 - phi)
        scale = np.abs(ref[..., 0]).max()
        deviation = np.abs(u - ref).max(axis=(0, 1)) / scale
        miss = np.abs(u - mirrored).max() / scale
        print(f"C. {label} beam reflected by Cox-Munk (n = {n_index}): max |VDISORT - closed form| / max I = {deviation}")
        print(f"   the mirrored relative azimuth misses by {miss:.3f}")
        assert deviation.max() < 1e-8, (label, deviation)
        assert miss > 1e4 * 1e-8, miss
        # Fresnel reflection makes no V from an unpolarized beam (M_VI = 0); the polarized beam exercises all four
        polarization = np.abs(ref).max(axis=(0, 1)) / scale
        print(f"   max |Q|, |U|, |V| / max I = {polarization[1:]}")
        assert polarization[1:3].min() > 1e-2 and (label == "unpolarized" or polarization[3] > 1e-2), polarization


if __name__ == "__main__":
    test_jones_rayleigh_matches_closed_form()
    test_polarized_beam_single_scattering()
    test_cox_munk_polarization_structure()
    test_beam_reflection_by_cox_munk()
    print("VDISORT polarized beam and surface tests passed")
