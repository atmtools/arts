#!/usr/bin/env python3

"""VDISORT with a polarized beam and a polarized surface.

The polarized comparisons with RT3 and RT4 never excite the Q, U and V
columns of the beam, because both have an unpolarized beam, and RT3 allows
only a Lambertian surface with a beam.  This test covers both through
pyarts3:

A. Single scattering of a polarized beam in an optically thin Rayleigh
   layer, with ARTS's Rayleigh phase-matrix Fourier modes
   (vdisort.scattering_optics of a GasScatterer).  VDISORT's diffuse operator
   is zero, so its formal solution is exactly the single-scattering integral,
   tabulated below from an independent evaluation with the dipole field
   E_out = E_in - k_out (k_out . E_in) projected on the meridional bases.
B. The polarization structure of VDISORT's raw polarized Cox-Munk
   reflection, M / M00, against tabulated values of the reflection of the
   fields by the specular facet (Fresnel amplitudes in the facet's s and p
   directions, projected on the meridional bases).  The scalar M00 (the facet
   slope distribution and shadowing) is DISORT's, validated against the
   Fortran DISORT problems.
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

    # ARTS's Rayleigh C^m and S^m of Z(mu_i, 0 <- -mu0, phi')
    species = arts.ArrayOfScatteringSpecies()
    species.add(arts.GasScatterer(arts.ConstantGasScattering(1e-30), arts.RayleighGasScattering(0.0)))
    atm = arts.AtmPoint()
    atm.pressure = 8e4
    atm.temperature = 260.0
    optics = vd.scattering_optics(species, atm, 89e9, mu, [-mu0], nfourier)
    cosine = np.asarray(optics.cosine)[:, None, :, 0]
    sine = np.asarray(optics.sine)[:, None, :, 0]

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

    # Upward at the top on stream 2 (mu = 0.237) and downward at the bottom on
    # stream 13 (mu = -0.763), at the azimuths phi_out
    n = nquad // 2
    up, down = 2, n + 5
    phi_out = np.array([0.0, 0.3, phi0, 1.9, np.pi, 4.4])

    def field(model, phi):
        u = np.asarray(model.u(tau=np.array([0.0, depth]), phi=phi))
        return np.concatenate((u[0, :, up], u[1, :, down]))

    ref = np.array([
        [4.96754443532156e-07, -5.28396544175680e-08, -4.22368131194178e-07, 1.47849145912888e-07],
        [4.68768270264573e-07, -7.87060830005066e-08, -3.28660402008651e-07, 1.87555666276326e-07],
        [4.00278493395939e-07, -5.75944213175187e-08, -1.66095830889410e-07, 2.07619788611763e-07],
        [3.03059610902693e-07, 8.30509129165863e-08, 2.80577308780512e-07, 4.55486824606329e-08],
        [5.78492366715659e-07, 2.88982687659360e-08, 3.99551092159033e-07, -2.40954513802322e-07],
        [5.13005852796478e-07, 2.38817457481129e-07, -5.12435970155258e-09, -2.62116365979781e-07],
        [1.67462078813691e-07, -3.47166480900049e-09, -7.36756627429650e-08, 8.68012905116077e-08],
        [1.85288199876436e-07, 1.50137353395687e-08, -8.37779116922669e-08, 9.50220435426892e-08],
        [1.96908969455559e-07, 5.45021792309593e-08, -7.93408612059474e-08, 9.91760765074343e-08],
        [1.71441864752456e-07, 1.03015164908675e-07, 7.65620969531610e-08, 6.56212214119098e-08],
        [1.13050884714876e-07, -5.78828589078150e-08, 9.64927349116139e-08, 6.30421257972797e-09],
        [5.22267701785750e-08, -3.30507992931378e-08, -4.03012092139639e-08, 1.92290797688661e-09],
    ])  # fmt: skip
    model = solve(beam)
    scale = np.abs(ref[:, 0]).max()
    deviation = np.abs(field(model, phi_out) - ref).max(axis=0) / scale
    print(f"A. polarized beam, single scattering: max |VDISORT - expected| / max I = {deviation}")
    assert deviation.max() < 1e-8, deviation

    # The comparison is not blind to the beam's polarization or the azimuth sense
    unpolarized = np.abs(field(solve(np.array([beam[0], 0, 0, 0])), phi_out) - ref).max() / scale
    mirrored = np.abs(field(model, np.mod(2 * phi0 - phi_out, 2 * np.pi)) - ref).max() / scale
    print(f"   an unpolarized beam misses by {unpolarized:.3f}, the mirrored azimuth by {mirrored:.3f}")
    assert unpolarized > 1e4 * 1e-8 and mirrored > 1e4 * 1e-8


def test_cox_munk_polarization_structure():
    """B. VDISORT's raw polarized Cox-Munk matrix M / M00 at one geometry per refractive index."""
    cases = {
        (1.34, 0.7, 0.5, 0.4): [
            [1.0, -0.83251026003135, 0.54322552353009, 0.0],
            [-0.75648400812752, 0.59875761966399, -0.47496322486588, 0.0],
            [0.6449025103993, -0.58855360289422, 0.28519572571752, 0.0],
            [0.0, 0.0, 0.0, -0.1087782033681],
        ],
        (1.5 + 0.3j, 0.35, 0.6, 5.7): [
            [1.0, -0.34307550922659, -0.89849786575612, 0.0],
            [-0.51052810861047, 0.1906592964045, 0.49540192622843, -0.23208984050515],
            [-0.81508298357774, 0.30148903141806, 0.79204359600663, 0.14536972270075],
            [0.0, 0.25584171251155, -0.09768863026451, 0.00165222783507],
        ],
        (7.0 + 3.0j, 0.9, 0.6, 0.9): [
            [1.0, -0.07413000690503, 0.06745388778762, 0.0],
            [-0.01232330780915, -0.57297917382974, -0.81238101640821, 0.0428780777427],
            [0.09946567742717, -0.81627157962926, 0.57751244399023, 0.00531238276314],
            [0.0, -0.02907830237946, -0.03195627156381, -0.99402613855827],
        ],
    }
    for (n, mu_out, mu_in, azimuth), expected in cases.items():
        raw = np.asarray(vd.cox_munk_reflection(mu_out, mu_in, azimuth, 7.0, complex(n), True))
        worst = np.abs(raw / raw[0, 0] - expected).max()
        print(f"B. Cox-Munk, n = {n}: max |raw / raw_00 - expected| = {worst:.1e}")
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
    test_polarized_beam_single_scattering()
    test_cox_munk_polarization_structure()
    test_beam_reflection_by_cox_munk()
    print("VDISORT polarized beam and surface tests passed")
