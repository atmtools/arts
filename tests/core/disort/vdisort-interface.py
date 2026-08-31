#!/usr/bin/env python3

"""Smoke tests for the low-level polarized VDISORT Python interface."""

import os

import numpy as np

from pyarts3 import arts


def _fresnel_amplitudes(mu, refractive_index):
    """Air-to-dielectric Fresnel amplitudes in the local (v, h) basis."""
    sin_incident = np.sqrt(1.0 - mu**2)
    sin_transmitted = sin_incident / refractive_index
    cos_transmitted = np.sqrt(1.0 - sin_transmitted**2)

    vertical = (refractive_index * mu - cos_transmitted) / (
        refractive_index * mu + cos_transmitted
    )
    horizontal = (mu - refractive_index * cos_transmitted) / (
        mu + refractive_index * cos_transmitted
    )
    return vertical, horizontal


def test_phase_helpers():
    cosine = np.zeros((1, 1, 2, 2, 4, 4))
    sine = np.zeros_like(cosine)
    cosine[0, 0, 0, 0, 0, 0] = 1.0

    combined = np.asarray(arts.vdisort.combine_phase_matrices(cosine, sine))
    assert combined.shape == (2, 1, 1, 2, 2, 4, 4)
    assert combined[arts.vdisort.cosine_mode, 0, 0, 0, 0, 0, 0] == 1.0

    beam_cosine = np.zeros((1, 1, 2, 4, 4))
    beam_sine = np.zeros_like(beam_cosine)
    combined_beam = np.asarray(
        arts.vdisort.combine_beam_phase_matrices(beam_cosine, beam_sine)
    )
    assert combined_beam.shape == (2, 1, 1, 2, 4, 4)


def test_bdrf():
    def cosine(mu_out, mu_in):
        return np.ones((4 * len(mu_out), 4 * len(mu_in)))

    bdrf = arts.VDisortBDRF(cosine)
    mu_out = np.array([0.2, 0.8])
    mu_in = np.array([0.4])
    assert np.asarray(bdrf(arts.vdisort.cosine_mode, mu_out, mu_in)).shape == (8, 4)
    np.testing.assert_array_equal(
        bdrf(arts.vdisort.sine_mode, mu_out, mu_in), np.zeros((8, 4))
    )


def test_absorbing_stokes_field():
    nquad = 2
    nfourier = 1
    depth = 0.7
    phase = np.zeros((2, nfourier, 1, nquad, nquad, 4, 4))
    b_pos = np.zeros((2, nfourier, nquad // 2, 4))
    b_neg = np.zeros_like(b_pos)

    top = np.array([2.0, -0.3, 0.2, -0.1])
    bottom = np.array([1.0, 0.25, -0.15, 0.05])
    b_neg[arts.vdisort.cosine_mode, 0, 0, :2] = top[:2]
    b_pos[arts.vdisort.cosine_mode, 0, 0, :2] = bottom[:2]
    b_neg[arts.vdisort.sine_mode, 0, 0, 2:] = top[2:]
    b_pos[arts.vdisort.sine_mode, 0, 0, 2:] = bottom[2:]

    model = arts.cppvdisort(
        tau_arr=np.array([depth]),
        omega_arr=np.array([0.0]),
        NQuad=nquad,
        phase_matrix=phase,
        mu0=0.5,
        beam_stokes=np.zeros(4),
        phi0=0.0,
        b_pos=b_pos,
        b_neg=b_neg,
    )

    tau = 0.3
    field = np.asarray(model.u(np.array([tau]), np.array([1.3])))
    assert field.shape == (1, 1, nquad, arts.vdisort.stokes_dimension)

    # DISORT uses one Gauss-Legendre rule per hemisphere (double Gauss),
    # whose one-point positive ordinate is the midpoint of [0, 1].
    mu = 0.5
    np.testing.assert_allclose(field[0, 0, 0], bottom * np.exp((tau - depth) / mu))
    np.testing.assert_allclose(field[0, 0, 1], top * np.exp(-tau / mu))

    compatible = np.asarray(model.pydisort_u(np.array([0.5, 0.2]), np.array([0.0])))
    assert compatible.shape == (nquad, 2, 1, arts.vdisort.stokes_dimension)
    assert np.asarray(model.flux(np.array([0.2, 0.5]))).shape == (3, 2)
    assert np.asarray(model.pydisort_flux_up(np.array([0.5, 0.2]))).shape == (2,)
    assert len(model.pydisort_flux_down(np.array([0.5, 0.2]))) == 2


def test_fresnel_brewster_angle():
    """Recover the air-to-glass Brewster minimum through VDISORT's boundary."""
    refractive_index = 1.5
    nquad = 24
    npositive = nquad // 2
    depth = 1.0e-8

    # VDISORT uses a Gauss-Legendre rule mapped from [-1, 1] to [0, 1]
    # independently in each hemisphere.
    nodes, weights = np.polynomial.legendre.leggauss(npositive)
    mu = 0.5 * (nodes + 1.0)
    weights = 0.5 * weights

    def specular_fresnel(mu_out, mu_in):
        """Mode-zero discrete specular BRDF, including quadrature normalization."""
        out = np.zeros((4 * len(mu_out), 4 * len(mu_in)))
        for i, outgoing in enumerate(mu_out):
            for j, incoming in enumerate(mu_in):
                if not np.isclose(outgoing, incoming, rtol=0.0, atol=1.0e-13):
                    continue
                rv, rh = _fresnel_amplitudes(incoming, refractive_index)
                reflection = np.asarray(arts.rtepack.fresnel_reflectance(rv, rh))
                out[4 * i : 4 * (i + 1), 4 * j : 4 * (j + 1)] = reflection / (
                    np.pi * weights[j] * incoming
                )
        return out

    phase = np.zeros((2, 1, 1, nquad, nquad, 4, 4))
    boundary_down = np.zeros((2, 1, npositive, 4))
    # Unit unpolarized illumination contains equal vertical (p) and horizontal
    # (s) intensities.  The two reflectances can therefore be recovered from
    # the reflected Stokes vector as I+Q and I-Q, respectively.
    boundary_down[arts.vdisort.cosine_mode, 0, :, 0] = 1.0

    model = arts.cppvdisort(
        tau_arr=np.array([depth]),
        omega_arr=np.array([0.0]),
        NQuad=nquad,
        phase_matrix=phase,
        mu0=0.5,
        beam_stokes=np.zeros(4),
        phi0=0.0,
        b_neg=boundary_down,
        BDRF_Fourier_modes=[arts.VDisortBDRF(specular_fresnel)],
    )

    stokes = np.asarray(model.u(np.array([depth]), np.array([0.0])))[0, 0, :npositive]
    rv, rh = _fresnel_amplitudes(mu, refractive_index)
    attenuation = np.exp(-depth / mu)
    expected_vertical = rv**2 * attenuation
    expected_horizontal = rh**2 * attenuation
    vertical = stokes[:, 0] + stokes[:, 1]
    horizontal = stokes[:, 0] - stokes[:, 1]

    # The full VDISORT solve must reproduce both Fresnel reflectances and must
    # not generate U or V for this coplanar case.
    np.testing.assert_allclose(vertical, expected_vertical, rtol=2.0e-8, atol=2.0e-12)
    np.testing.assert_allclose(
        horizontal, expected_horizontal, rtol=2.0e-8, atol=2.0e-12
    )
    np.testing.assert_allclose(stokes[:, 2:], 0.0, atol=2.0e-12)

    numerical_index = np.argmin(vertical)
    numerical_angle = np.arccos(mu[numerical_index])
    brewster_angle = np.arctan(refractive_index)

    assert numerical_index == np.argmin(np.abs(np.arccos(mu) - brewster_angle))
    assert abs(numerical_angle - brewster_angle) < np.deg2rad(0.6)
    assert vertical[numerical_index] < 5.0e-5
    # The orthogonal Fresnel component remains finite at the Brewster angle.
    assert horizontal[numerical_index] > 0.1

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        angle = np.linspace(0.0, 0.5 * np.pi - 1.0e-6, 500)
        dense_rv, dense_rh = _fresnel_amplitudes(np.cos(angle), refractive_index)
        quadrature_angle = np.arccos(mu)

        _, (reflectance_plot, error_plot) = plt.subplots(
            1, 2, figsize=(11, 4), constrained_layout=True
        )
        reflectance_plot.plot(np.rad2deg(angle), dense_rv**2, label="Fresnel p")
        reflectance_plot.plot(np.rad2deg(angle), dense_rh**2, label="Fresnel s")
        reflectance_plot.plot(
            np.rad2deg(quadrature_angle),
            vertical,
            "o",
            fillstyle="none",
            label="VDISORT p",
        )
        reflectance_plot.plot(
            np.rad2deg(quadrature_angle),
            horizontal,
            "x",
            label="VDISORT s",
        )
        reflectance_plot.axvline(
            np.rad2deg(brewster_angle),
            color="black",
            linestyle="--",
            label="Brewster angle",
        )
        reflectance_plot.set(
            xlabel="Incidence angle [degree]",
            ylabel="Reflectance",
            title="Air-to-glass Fresnel reflection",
            xlim=(0.0, 90.0),
            ylim=(-0.02, 1.02),
        )
        reflectance_plot.legend()
        reflectance_plot.grid(True, alpha=0.25)

        error_plot.semilogy(
            np.rad2deg(quadrature_angle),
            np.maximum(np.abs(vertical - expected_vertical), 1.0e-18),
            "o-",
            label="p error",
        )
        error_plot.semilogy(
            np.rad2deg(quadrature_angle),
            np.maximum(np.abs(horizontal - expected_horizontal), 1.0e-18),
            "x-",
            label="s error",
        )
        error_plot.set(
            xlabel="Incidence angle [degree]",
            ylabel="Absolute error",
            title="VDISORT minus analytical Fresnel",
            xlim=(0.0, 90.0),
        )
        error_plot.legend()
        error_plot.grid(True, which="both", alpha=0.25)
        plt.show()


if __name__ == "__main__":
    test_phase_helpers()
    test_bdrf()
    test_absorbing_stokes_field()
    test_fresnel_brewster_angle()
