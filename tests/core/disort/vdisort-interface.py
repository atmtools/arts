#!/usr/bin/env python3

"""Smoke tests for the low-level polarized VDISORT Python interface."""

import numpy as np

from pyarts3 import arts


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


if __name__ == "__main__":
    test_phase_helpers()
    test_bdrf()
    test_absorbing_stokes_field()
