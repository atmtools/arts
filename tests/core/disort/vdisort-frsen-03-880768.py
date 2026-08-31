#!/usr/bin/env python3

"""Self-contained VDISORT tests derived from Lin et al. (2022), 880768.

The source is ``doi: 10.3389/frsen.2022.880768``
"""

import os

import numpy as np

from pyarts3 import arts

COSINE = arts.vdisort.cosine_mode
SINE = arts.vdisort.sine_mode


def _double_gauss_mu(nquad):
    nodes, _ = np.polynomial.legendre.leggauss(nquad // 2)
    positive = 0.5 * (nodes + 1.0)
    return np.concatenate((positive, -positive))


def _empty_phase(nfourier, nlayers, nquad):
    return np.zeros((2, nfourier, nlayers, nquad, nquad, 4, 4))


def _missing_figure(figure, reason, **missing_data):
    """Report an intentionally empty paper-figure comparison.

    Zero-length arrays mark data that must be replaced by published or
    author-supplied values.  They must never be interpreted as physical zero.
    """
    names = ", ".join(missing_data)
    if any(np.asarray(value).size for value in missing_data.values()):
        raise ValueError(f"Figure {figure} placeholder data must remain empty")
    print(f"Figure {figure} comparison is not implemented: {reason}")
    print(f"Required data: {names}")


def test_layered_single_scattering_recurrence():
    """Check the two-layer specialization of paper Eqs. (136)--(141)."""
    interfaces = np.array([0.17, 0.53])
    omega = np.array([0.35, 0.82])
    phase_value = np.array([0.7, 1.6])
    mu0 = 0.73
    user_mu = np.linspace(0.05, 0.95, 61)
    nquad = 4

    phase = _empty_phase(1, 2, nquad)
    beam_phase = np.zeros((2, 1, 2, nquad, 4, 4))
    beam_phase[COSINE, 0, :, :, 0, 0] = phase_value[:, None]
    model = arts.cppvdisort(
        tau_arr=interfaces,
        omega_arr=omega,
        NQuad=nquad,
        phase_matrix=phase,
        mu0=mu0,
        beam_stokes=np.array([1.0, 0.0, 0.0, 0.0]),
        phi0=0.0,
        beam_phase_matrix=beam_phase,
    )

    user_phase = np.zeros((2, 1, 2, len(user_mu), nquad, 4, 4))
    user_beam = np.zeros((2, 1, 2, len(user_mu), 4, 4))
    user_beam[COSINE, 0, :, :, 0, 0] = phase_value[:, None]
    result = np.asarray(
        model.u_user(
            tau=np.array([0.0]),
            phi=np.array([0.0]),
            mu=user_mu,
            phase_matrix=user_phase,
            beam_phase_matrix=user_beam,
        )
    )[0, 0]

    expected = np.zeros_like(user_mu)
    tops = np.concatenate(([0.0], interfaces[:-1]))
    for top, bottom, layer_omega, layer_phase in zip(
        tops, interfaces, omega, phase_value, strict=True
    ):
        extinction = 1.0 / mu0 + 1.0 / user_mu
        expected += (
            layer_omega
            * layer_phase
            * (np.exp(-top * extinction) - np.exp(-bottom * extinction))
            / (4.0 * np.pi * user_mu * extinction)
        )

    np.testing.assert_allclose(result[:, 0], expected, rtol=3.0e-9, atol=2.0e-14)
    np.testing.assert_allclose(result[:, 1:], 0.0, atol=2.0e-14)

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        zenith = np.rad2deg(np.arccos(user_mu))
        figure, (radiance_plot, error_plot) = plt.subplots(
            1, 2, figsize=(11, 4), constrained_layout=True
        )
        radiance_plot.plot(zenith, result[:, 0], label="VDISORT")
        radiance_plot.plot(zenith, expected, "--", label="Eqs. 136--141")
        radiance_plot.set(
            xlabel="Upward zenith angle [degree]",
            ylabel="Radiance I",
            title="Two-layer single scattering",
        )
        radiance_plot.grid(True, alpha=0.25)
        radiance_plot.legend()
        error_plot.semilogy(
            zenith, np.maximum(np.abs(result[:, 0] - expected), 1.0e-18)
        )
        error_plot.set(
            xlabel="Upward zenith angle [degree]",
            ylabel="Absolute error",
            title="VDISORT minus analytical recurrence",
        )
        error_plot.grid(True, which="both", alpha=0.25)
        plt.show()


def test_polarized_lambertian_depolarizes():
    """Verify Eq. (146), transformed from [Iv,Ih,U,V] to [I,Q,U,V]."""
    nquad = 8
    npositive = nquad // 2
    rho = 0.37
    depth = 1.0e-9
    incident = np.array([1.2, -0.31, 0.23, -0.17])

    boundary_down = np.zeros((2, 1, npositive, 4))
    boundary_down[COSINE, 0, :, :2] = incident[:2]
    boundary_down[SINE, 0, :, 2:] = incident[2:]

    def lambertian(mu_out, mu_in):
        blocks = np.zeros((len(mu_out), 4, len(mu_in), 4))
        blocks[:, 0, :, 0] = rho / np.pi
        return blocks.reshape(4 * len(mu_out), 4 * len(mu_in))

    model = arts.cppvdisort(
        tau_arr=np.array([depth]),
        omega_arr=np.array([0.0]),
        NQuad=nquad,
        phase_matrix=_empty_phase(1, 1, nquad),
        mu0=0.5,
        beam_stokes=np.zeros(4),
        phi0=0.0,
        b_neg=boundary_down,
        BDRF_Fourier_modes=[arts.VDisortBDRF(lambertian)],
    )

    field = np.asarray(model.u(np.array([depth]), np.array([0.0, 1.1, 2.7])))[0]
    upward = field[:, :npositive]
    # Integral_0^1 mu dmu = 1/2; the transparent-layer attenuation tends to one.
    expected_i = 0.5 * rho * incident[0]
    np.testing.assert_allclose(upward[..., 0], expected_i, rtol=2.0e-8)
    np.testing.assert_allclose(upward[..., 1:], 0.0, atol=2.0e-12)
    np.testing.assert_allclose(
        upward[1:], np.broadcast_to(upward[0], upward[1:].shape), atol=2.0e-12
    )

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        figure, (stokes_plot, angular_plot) = plt.subplots(
            1, 2, figsize=(11, 4), constrained_layout=True
        )
        components = np.arange(4)
        width = 0.35
        stokes_plot.bar(components - width / 2, incident, width, label="Incident")
        stokes_plot.bar(
            components + width / 2,
            upward[0, 0],
            width,
            label="Lambertian reflection",
        )
        stokes_plot.set(
            xticks=components,
            xticklabels=list("IQUV"),
            ylabel="Stokes component",
            title="Equation 146 depolarization",
        )
        stokes_plot.axhline(0.0, color="black", linewidth=0.8)
        stokes_plot.legend()

        for component, label in enumerate("IQUV"):
            angular_plot.plot(
                [0.0, 63.0, 155.0], upward[:, 0, component], "o-", label=label
            )
        angular_plot.axhline(
            expected_i, color="black", linestyle="--", label="Analytical I"
        )
        angular_plot.set(
            xlabel="Azimuth [degree]",
            ylabel="Reflected radiance",
            title="Azimuth-independent Lambertian field",
        )
        angular_plot.grid(True, alpha=0.25)
        angular_plot.legend()
        plt.show()


def _decoupled_model(v_block):
    nquad = 2
    phase = _empty_phase(1, 1, nquad)
    for outgoing in range(nquad):
        for incoming in range(nquad):
            phase[COSINE, 0, 0, outgoing, incoming, :2, :2] = np.array(
                [[0.62, 0.11], [0.11, 0.47]]
            )
            phase[SINE, 0, 0, outgoing, incoming, 2, 2] = 0.38
            phase[SINE, 0, 0, outgoing, incoming, 3, 3] = v_block
    down = np.zeros((2, 1, 1, 4))
    down[COSINE, 0, 0, :2] = [1.0, -0.2]
    down[SINE, 0, 0, 2] = 0.16
    return arts.cppvdisort(
        tau_arr=np.array([0.6]),
        omega_arr=np.array([0.7]),
        NQuad=nquad,
        phase_matrix=phase,
        mu0=0.5,
        beam_stokes=np.zeros(4),
        phi0=0.0,
        b_neg=down,
    )


def test_b2_zero_decouples_circular_polarization():
    """Exercise the real 3x3 limit obtained when the paper's b2 is zero."""
    first = _decoupled_model(0.15)
    second = _decoupled_model(0.71)
    assert not first.has_complex_eigensolutions()
    assert not second.has_complex_eigensolutions()

    tau = np.linspace(0.0, 0.6, 41)
    phi = np.array([0.0, 0.8])
    first_field = np.asarray(first.u(tau, phi))
    second_field = np.asarray(second.u(tau, phi))
    np.testing.assert_allclose(first_field[..., :3], second_field[..., :3], atol=2e-12)
    np.testing.assert_allclose(first_field[..., 3], 0.0, atol=2e-12)
    np.testing.assert_allclose(second_field[..., 3], 0.0, atol=2e-12)

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        figure, (field_plot, difference_plot) = plt.subplots(
            1, 2, figsize=(11, 4), constrained_layout=True
        )
        for component, label in enumerate("IQU"):
            field_plot.plot(tau, first_field[:, 0, 0, component], label=label)
        field_plot.set(
            xlabel="Optical depth",
            ylabel="Upward Stokes component",
            title=r"Real $3\times3$ limit ($b_2=0$)",
        )
        field_plot.grid(True, alpha=0.25)
        field_plot.legend()
        difference_plot.semilogy(
            tau,
            np.maximum(
                np.max(
                    np.abs(first_field[..., :3] - second_field[..., :3]), axis=(1, 2, 3)
                ),
                1.0e-18,
            ),
            label="max |delta(I,Q,U)|",
        )
        difference_plot.semilogy(
            tau,
            np.maximum(np.max(np.abs(first_field[..., 3]), axis=(1, 2)), 1.0e-18),
            label="max |V|",
        )
        difference_plot.set(
            xlabel="Optical depth",
            ylabel="Absolute invariant residual",
            title="Independence from the decoupled V block",
        )
        difference_plot.grid(True, which="both", alpha=0.25)
        difference_plot.legend()
        plt.show()


def _local_basis(mu, phi):
    direction = np.array(
        [
            np.sqrt(max(0.0, 1.0 - mu * mu)) * np.cos(phi),
            np.sqrt(max(0.0, 1.0 - mu * mu)) * np.sin(phi),
            mu,
        ]
    )
    horizontal = np.array([np.sin(phi), -np.cos(phi), 0.0])
    vertical = np.cross(horizontal, direction)
    return direction, vertical, horizontal


def _mueller_from_real_jones(jones):
    stokes_coherencies = (
        np.array([[0.5, 0.0], [0.0, 0.5]], dtype=complex),
        np.array([[0.5, 0.0], [0.0, -0.5]], dtype=complex),
        np.array([[0.0, 0.5], [0.5, 0.0]], dtype=complex),
        np.array([[0.0, -0.5j], [0.5j, 0.0]], dtype=complex),
    )
    out = np.empty((4, 4))
    for column, coherency in enumerate(stokes_coherencies):
        scattered = jones @ coherency @ jones.T
        out[:, column] = (
            np.trace(scattered).real,
            (scattered[0, 0] - scattered[1, 1]).real,
            (scattered[0, 1] + scattered[1, 0]).real,
            (-1j * (scattered[0, 1] - scattered[1, 0])).real,
        )
    return 1.5 * out


def test_figure_3_cox_munk_reflection():
    """Reproduce paper Figure 3 from its analytic 5 m/s ocean BPrDF."""
    nquad = 2
    nfourier = 24
    mu0 = np.cos(np.deg2rad(30.0))
    phase = _empty_phase(nfourier, 1, nquad)
    modes = arts.vdisort.cox_munk_fourier_modes(
        wind_speed=5.0,
        refractive_index=1.34,
        shadowing=True,
        number_of_modes=nfourier,
        azimuth_quadrature_points=64,
    )
    model = arts.cppvdisort(
        tau_arr=np.array([1.0e-10]),
        omega_arr=np.array([0.0]),
        NQuad=nquad,
        NFourier=nfourier,
        phase_matrix=phase,
        mu0=mu0,
        beam_stokes=np.array([1.0, 0.0, 0.0, 0.0]),
        phi0=0.0,
        BDRF_Fourier_modes=modes,
    )

    polar_angle = np.deg2rad(np.linspace(2.0, 82.0, 41))
    mu = np.cos(polar_angle)
    azimuth = np.deg2rad(np.array([0.0, 45.0, 90.0, 135.0, 180.0]))
    user_phase = np.zeros((2, nfourier, 1, len(mu), nquad, 4, 4))
    user_beam_phase = np.zeros((2, nfourier, 1, len(mu), 4, 4))
    reflected = np.asarray(
        model.u_user(
            tau=np.array([1.0e-10]),
            phi=azimuth,
            mu=mu,
            phase_matrix=user_phase,
            beam_phase_matrix=user_beam_phase,
        )
    )[0]
    analytic = np.array(
        [
            [
                mu0
                * (
                    np.asarray(
                        arts.vdisort.cox_munk_reflection(
                            out_mu, mu0, phi, 5.0, 1.34, True
                        )
                    )
                    @ np.array([1, 0, 0, 0])
                )
                for out_mu in mu
            ]
            for phi in azimuth
        ]
    )

    # A finite Fourier truncation is least accurate in the narrow specular
    # peak.  The integral shape and all three nonzero Stokes components must
    # nevertheless track the analytic BPrDF.
    scale = np.maximum(np.abs(analytic[..., :3]), 2.0e-4)
    np.testing.assert_allclose(reflected[..., 3], 0.0, atol=2.0e-10)
    maximum_relative_error = np.max(
        np.abs(reflected[..., :3] - analytic[..., :3]) / scale
    )
    assert maximum_relative_error < 0.12, maximum_relative_error

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots(3, 5, figsize=(14, 8), constrained_layout=True)
        for column, angle in enumerate(np.rad2deg(azimuth)):
            for row, (component, label) in enumerate(zip(range(3), "IQU", strict=True)):
                axes[row, column].plot(
                    np.rad2deg(polar_angle),
                    analytic[column, :, component],
                    "r-",
                    label="Analytic",
                )
                axes[row, column].plot(
                    np.rad2deg(polar_angle),
                    reflected[column, :, component],
                    "ko",
                    markersize=3,
                    fillstyle="none",
                    label="VDISORT",
                )
                axes[row, column].set_title(rf"$\phi={angle:.0f}^\circ$")
                axes[row, column].grid(True, alpha=0.2)
                if column == 0:
                    axes[row, column].set_ylabel(label)
                if row == 2:
                    axes[row, column].set_xlabel("Polar angle [degree]")
        axes[0, 0].legend()
        figure.suptitle("Reproduction of Lin et al. Figure 3: 5 m/s Cox-Munk surface")
        plt.show()


def _rayleigh_lab_matrix(mu_out, phi_out, mu_in, phi_in):
    _, out_v, out_h = _local_basis(mu_out, phi_out)
    _, in_v, in_h = _local_basis(mu_in, phi_in)
    jones = np.array(
        [
            [np.dot(out_v, in_v), np.dot(out_v, in_h)],
            [np.dot(out_h, in_v), np.dot(out_h, in_h)],
        ]
    )
    return _mueller_from_real_jones(jones)


def _rayleigh_fourier(mu_out, mu_in, nfourier, samples=128):
    """Numerically form the ordinary Fourier matrices, then combine them."""
    phi = 2.0 * np.pi * np.arange(samples) / samples
    cosine = np.zeros((nfourier, 1, len(mu_out), len(mu_in), 4, 4))
    sine = np.zeros_like(cosine)
    for outgoing, out_mu in enumerate(mu_out):
        for incoming, in_mu in enumerate(mu_in):
            matrices = np.array(
                [_rayleigh_lab_matrix(out_mu, 0.0, in_mu, angle) for angle in phi]
            )
            for mode in range(nfourier):
                epsilon = 1.0 if mode == 0 else 2.0
                cosine[mode, 0, outgoing, incoming] = epsilon * np.mean(
                    matrices * np.cos(mode * phi)[:, None, None], axis=0
                )
                sine[mode, 0, outgoing, incoming] = epsilon * np.mean(
                    matrices * np.sin(mode * phi)[:, None, None], axis=0
                )
    return np.asarray(arts.vdisort.combine_phase_matrices(cosine, sine))


def _rayleigh_beam_fourier(mu_out, mu0, nfourier, samples=128):
    """Project an unpolarized beam using VDISORT's combined-vector layout."""
    phi = 2.0 * np.pi * np.arange(samples) / samples
    out = np.zeros((2, nfourier, 1, len(mu_out), 4, 4))
    for outgoing, out_mu in enumerate(mu_out):
        column = np.array(
            [_rayleigh_lab_matrix(out_mu, angle, -mu0, 0.0)[:, 0] for angle in phi]
        )
        delta = -phi
        for mode in range(nfourier):
            epsilon = 1.0 if mode == 0 else 2.0
            cosine = epsilon * np.mean(column * np.cos(mode * delta)[:, None], axis=0)
            sine = epsilon * np.mean(column * np.sin(mode * delta)[:, None], axis=0)
            out[COSINE, mode, 0, outgoing, :2, 0] = cosine[:2] / epsilon
            out[SINE, mode, 0, outgoing, :2, 0] = sine[:2] / epsilon
            out[SINE, mode, 0, outgoing, 2:, 0] = cosine[2:] / epsilon
            out[COSINE, mode, 0, outgoing, 2:, 0] = sine[2:] / epsilon
    return out


def test_multiple_scattering_rayleigh_principal_plane():
    """Check the principal-plane Rayleigh symmetries used in paper section 6."""
    nquad = 8
    nfourier = 3
    mu0 = 0.6
    quadrature_mu = _double_gauss_mu(nquad)
    phase = _rayleigh_fourier(quadrature_mu, quadrature_mu, nfourier)
    beam = _rayleigh_beam_fourier(quadrature_mu, mu0, nfourier)

    model = arts.cppvdisort(
        tau_arr=np.array([0.32]),
        omega_arr=np.array([0.95]),
        NQuad=nquad,
        NFourier=nfourier,
        phase_matrix=phase,
        mu0=mu0,
        beam_stokes=np.array([1.0, 0.0, 0.0, 0.0]),
        phi0=0.0,
        beam_phase_matrix=beam,
    )
    user_mu = np.linspace(0.05, 0.95, 37)
    user_phase = _rayleigh_fourier(user_mu, quadrature_mu, nfourier)
    user_beam = _rayleigh_beam_fourier(user_mu, mu0, nfourier)
    field = np.asarray(
        model.u_user(
            tau=np.array([0.0]),
            phi=np.array([0.0, np.pi]),
            mu=user_mu,
            phase_matrix=user_phase,
            beam_phase_matrix=user_beam,
        )
    )[0]

    assert np.all(np.isfinite(field))
    assert np.max(np.abs(field[..., 1])) > 1.0e-5
    np.testing.assert_allclose(field[..., 2:], 0.0, atol=2.0e-11)

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        zenith = np.rad2deg(np.arccos(user_mu))
        figure, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
        titles = ("Forward principal plane", "Backward principal plane")
        for plane, title in enumerate(titles):
            axes[0, plane].plot(zenith, field[plane, :, 0], label="I")
            axes[0, plane].plot(zenith, field[plane, :, 1], label="Q")
            axes[0, plane].set(
                ylabel="Top-of-atmosphere radiance",
                title=title,
            )
            axes[0, plane].grid(True, alpha=0.25)
            axes[0, plane].legend()
            axes[1, plane].plot(
                zenith, field[plane, :, 1] / field[plane, :, 0], label="Q/I"
            )
            axes[1, plane].plot(
                zenith, field[plane, :, 2] / field[plane, :, 0], label="U/I"
            )
            axes[1, plane].plot(
                zenith, field[plane, :, 3] / field[plane, :, 0], label="V/I"
            )
            axes[1, plane].set(
                xlabel="Upward zenith angle [degree]",
                ylabel="Normalized polarization",
            )
            axes[1, plane].grid(True, alpha=0.25)
            axes[1, plane].legend()
        figure.suptitle("Rayleigh multiple scattering: paper principal-plane symmetry")
        plt.show()


def test_figure_2_missing_garcia_siewert_data():
    """Placeholder for Figure 2's Garcia-Siewert benchmark comparison."""
    _missing_figure(
        2,
        "the article identifies the benchmark but does not tabulate its inputs or curves",
        polarized_phase_expansion=np.zeros((0, 4, 4)),
        benchmark_angles=np.zeros(0),
        benchmark_stokes=np.zeros((0, 4)),
    )


def test_figure_4_missing_cirrus_benchmark_data():
    """Placeholder for the Figure 4 cirrus/VLBLE comparison."""
    _missing_figure(
        4,
        "the b2=0 invariant is tested separately, but the cirrus expansion and VLBLE curve are external",
        cirrus_phase_expansion=np.zeros((0, 4, 4)),
        vlble_angles=np.zeros(0),
        vlble_stokes=np.zeros((0, 4)),
    )


def test_figure_5_missing_ice_cloud_phase_data():
    """Placeholder for Figure 5's ice-cloud phase-matrix elements."""
    _missing_figure(
        5,
        "the required Kokhanovsky et al. ice-cloud phase matrix is not supplied by the article",
        ice_cloud_scattering_angles=np.zeros(0),
        ice_cloud_mueller_matrices=np.zeros((0, 4, 4)),
    )


def _missing_polarized_aerosol_figure(figure):
    _missing_figure(
        figure,
        "bulk aerosol properties are stated, but neither the Mie expansion nor the plotted reference values are supplied",
        aerosol_phase_expansion=np.zeros((0, 4, 4)),
        reference_angles=np.zeros(0),
        reference_stokes=np.zeros((0, 4)),
    )


def test_figure_6_missing_polarized_aerosol_data():
    """Placeholder for the Figure 6 aerosol comparison."""
    _missing_polarized_aerosol_figure(6)


def test_figure_7_missing_polarized_aerosol_data():
    """Placeholder for the Figure 7 aerosol comparison."""
    _missing_polarized_aerosol_figure(7)


def test_figure_8_missing_polarized_aerosol_data():
    """Placeholder for the Figure 8 aerosol comparison."""
    _missing_polarized_aerosol_figure(8)


def _missing_polarized_cloud_figure(figure):
    _missing_figure(
        figure,
        "bulk cloud properties are stated, but neither the Mie expansion nor the plotted reference values are supplied",
        cloud_phase_expansion=np.zeros((0, 4, 4)),
        reference_angles=np.zeros(0),
        reference_stokes=np.zeros((0, 4)),
    )


def test_figure_9_missing_polarized_cloud_data():
    """Placeholder for the Figure 9 cloud comparison."""
    _missing_polarized_cloud_figure(9)


def test_figure_10_missing_polarized_cloud_data():
    """Placeholder for the Figure 10 cloud comparison."""
    _missing_polarized_cloud_figure(10)


def test_figure_11_missing_coupled_atmosphere_surface_reference():
    """Placeholder for Figure 11's atmosphere-over-ocean comparison."""
    _missing_figure(
        11,
        "the geometry is stated and the ocean is tested in Figure 3, but the molecular depolarization convention and reference curves are unavailable",
        molecular_phase_expansion=np.zeros((0, 4, 4)),
        reference_angles=np.zeros(0),
        reference_stokes=np.zeros((0, 4)),
    )


def _missing_ssc_aerosol_figure(figure):
    _missing_figure(
        figure,
        "the aerosol expansion, paper-specific single-scattering-correction inputs, and benchmark curves are unavailable",
        aerosol_phase_expansion=np.zeros((0, 4, 4)),
        untruncated_phase_matrix=np.zeros((0, 4, 4)),
        ssc_reference_angles=np.zeros(0),
        ssc_reference_stokes=np.zeros((0, 4)),
    )


def test_figure_12_missing_ssc_aerosol_data():
    """Placeholder for the Figure 12 SSC comparison."""
    _missing_ssc_aerosol_figure(12)


def test_figure_13_missing_ssc_aerosol_data():
    """Placeholder for the Figure 13 SSC comparison."""
    _missing_ssc_aerosol_figure(13)


def test_figure_14_missing_ssc_aerosol_data():
    """Placeholder for the Figure 14 SSC comparison."""
    _missing_ssc_aerosol_figure(14)


def test_figure_15_missing_l13_phase_data():
    """Placeholder for Figure 15's L=13 benchmark comparison."""
    _missing_figure(
        15,
        "the scalar inputs are stated, but the external L=13 phase matrix and Schulz et al. reference arrays are unavailable",
        l13_phase_expansion=np.zeros((0, 4, 4)),
        schulz_angles=np.zeros(0),
        schulz_stokes=np.zeros((0, 4)),
    )


if __name__ == "__main__":
    test_layered_single_scattering_recurrence()
    test_polarized_lambertian_depolarizes()
    test_b2_zero_decouples_circular_polarization()
    test_figure_3_cox_munk_reflection()
    test_multiple_scattering_rayleigh_principal_plane()
    test_figure_2_missing_garcia_siewert_data()
    test_figure_4_missing_cirrus_benchmark_data()
    test_figure_5_missing_ice_cloud_phase_data()
    test_figure_6_missing_polarized_aerosol_data()
    test_figure_7_missing_polarized_aerosol_data()
    test_figure_8_missing_polarized_aerosol_data()
    test_figure_9_missing_polarized_cloud_data()
    test_figure_10_missing_polarized_cloud_data()
    test_figure_11_missing_coupled_atmosphere_surface_reference()
    test_figure_12_missing_ssc_aerosol_data()
    test_figure_13_missing_ssc_aerosol_data()
    test_figure_14_missing_ssc_aerosol_data()
    test_figure_15_missing_l13_phase_data()
