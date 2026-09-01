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


def _rayleigh_phase_matrix(cos_scattering_angle):
    """Normalized Rayleigh Mueller matrix in the scattering-plane basis."""
    cosine = np.asarray(cos_scattering_angle)
    phase = np.zeros(cosine.shape + (4, 4))
    phase[..., 0, 0] = 0.75 * (1.0 + cosine**2)
    phase[..., 0, 1] = -0.75 * (1.0 - cosine**2)
    phase[..., 1, 0] = phase[..., 0, 1]
    phase[..., 1, 1] = phase[..., 0, 0]
    phase[..., 2, 2] = 1.5 * cosine
    phase[..., 3, 3] = 1.5 * cosine
    return phase


def _sun_halo_profile(scattering_angle, center=np.deg2rad(22.0)):
    """Smooth axisymmetric phase profile with a narrow 22-degree halo."""
    width = np.deg2rad(1.25)
    return 1.0 + 80.0 * np.exp(-0.5 * ((scattering_angle - center) / width) ** 2)


def test_combined_surface_models():
    """Compose polarized surface operators after their Fourier projection."""
    mu = np.array([0.2, 0.7])
    lambertian_albedo = 0.8
    fresnel_fraction = 0.02
    number_of_modes = 3

    lambertian = arts.vdisort.lambertian_fourier_modes(
        lambertian_albedo, number_of_modes
    )
    fresnel = arts.vdisort.fresnel_fourier_modes(1.5, number_of_modes)
    generic = arts.vdisort.combine_fourier_modes(
        fresnel, fresnel_fraction, lambertian, 1.0 - fresnel_fraction
    )
    named = arts.vdisort.fresnel_lambertian_fourier_modes(
        fresnel_fraction, lambertian_albedo, 1.5, number_of_modes
    )

    for mode in range(number_of_modes):
        for alpha in (arts.vdisort.cosine_mode, arts.vdisort.sine_mode):
            expected = fresnel_fraction * np.asarray(fresnel[mode](alpha, mu, mu))
            expected += (1.0 - fresnel_fraction) * np.asarray(
                lambertian[mode](alpha, mu, mu)
            )
            np.testing.assert_allclose(generic[mode](alpha, mu, mu), expected)
            np.testing.assert_allclose(named[mode](alpha, mu, mu), expected)

    first_lambertian = np.asarray(
        lambertian[0](arts.vdisort.cosine_mode, mu, mu)
    ).reshape(2, 4, 2, 4)
    np.testing.assert_allclose(
        first_lambertian[:, 0, :, 0], 2.0 * lambertian_albedo / np.pi
    )
    first_lambertian[:, 0, :, 0] = 0.0
    np.testing.assert_allclose(first_lambertian, 0.0)

    rough = arts.vdisort.cox_munk_lambertian_fourier_modes(
        0.35, lambertian_albedo, 5.0, 1.34, True, 2, 24
    )
    assert len(rough) == 2
    assert np.isfinite(rough[0](arts.vdisort.cosine_mode, mu, mu)).all()

    scalar = arts.disort.cox_munk_lambertian_fourier_modes(
        0.35, lambertian_albedo, 5.0, 1.34, True, 2, 24
    )
    assert len(scalar) == 2
    assert np.isfinite(scalar[0](mu, mu)).all()


def _rayleigh_stokes_column(mu_out, phi_out, mu0, phi0):
    """Rayleigh first Mueller column in each outgoing local-meridian frame."""
    mu_out, phi_out = np.broadcast_arrays(mu_out, phi_out)
    horizontal = np.sqrt(1.0 - mu_out**2)
    outgoing = np.stack(
        (
            horizontal * np.cos(phi_out),
            horizontal * np.sin(phi_out),
            mu_out,
        ),
        axis=-1,
    )
    incident = np.array(
        [
            np.sqrt(1.0 - mu0**2) * np.cos(phi0),
            np.sqrt(1.0 - mu0**2) * np.sin(phi0),
            -mu0,
        ]
    )
    cos_scattering = np.clip(outgoing @ incident, -1.0, 1.0)
    polarized = 0.75 * (1.0 - cos_scattering**2)

    # The local horizontal axis is k_out x z and the local vertical axis is
    # horizontal x k_out.  Rayleigh polarization is perpendicular to the
    # scattering plane spanned by the incident and outgoing rays.
    local_horizontal = np.stack(
        (np.sin(phi_out), -np.cos(phi_out), np.zeros_like(phi_out)), axis=-1
    )
    local_vertical = np.cross(local_horizontal, outgoing)
    perpendicular = np.cross(incident, outgoing)
    perpendicular /= np.linalg.norm(perpendicular, axis=-1)[..., None]
    vertical_component = np.sum(perpendicular * local_vertical, axis=-1)
    horizontal_component = np.sum(perpendicular * local_horizontal, axis=-1)

    stokes = np.zeros(mu_out.shape + (4,))
    stokes[..., 0] = 0.75 * (1.0 + cos_scattering**2)
    stokes[..., 1] = polarized * (vertical_component**2 - horizontal_component**2)
    stokes[..., 2] = 2.0 * polarized * vertical_component * horizontal_component
    return stokes


def _rayleigh_beam_fourier(mu_out, mu0, phi0, nfourier):
    """Project the Rayleigh beam source into VDISORT's combined Fourier modes."""
    sample_phi = np.linspace(0.0, 2.0 * np.pi, 2048, endpoint=False)
    mu_grid, phi_grid = np.meshgrid(mu_out, sample_phi, indexing="ij")
    stokes = _rayleigh_stokes_column(mu_grid, phi_grid, mu0, phi0)
    delta_phi = phi0 - sample_phi
    combined = np.zeros((2, nfourier, 1, len(mu_out), 4, 4))

    for mode in range(nfourier):
        epsilon = 1.0 if mode == 0 else 2.0
        cosine = epsilon * np.mean(
            stokes * np.cos(mode * delta_phi)[None, :, None], axis=1
        )
        sine = epsilon * np.mean(
            stokes * np.sin(mode * delta_phi)[None, :, None], axis=1
        )
        # The beam source multiplies nonzero Fourier modes by epsilon.  Divide
        # the standard Fourier coefficients by the same factor here.  I,Q use
        # the ordinary combined ordering; U,V swap cosine and sine systems.
        combined[arts.vdisort.cosine_mode, mode, 0, :, :2, 0] = cosine[:, :2] / epsilon
        combined[arts.vdisort.sine_mode, mode, 0, :, :2, 0] = sine[:, :2] / epsilon
        combined[arts.vdisort.sine_mode, mode, 0, :, 2:, 0] = cosine[:, 2:] / epsilon
        combined[arts.vdisort.cosine_mode, mode, 0, :, 2:, 0] = sine[:, 2:] / epsilon
    return combined


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
    assert np.asarray(model.flux(np.array([0.2, 0.5]))).shape == (4, 2)


def test_fresnel_brewster_angle():
    """Recover the air-to-glass Brewster minimum through VDISORT's boundary."""
    refractive_index = 1.5
    nquad = 24
    npositive = nquad // 2
    depth = 1.0e-8

    # VDISORT uses a Gauss-Legendre rule mapped from [-1, 1] to [0, 1]
    # independently in each hemisphere.
    nodes, _ = np.polynomial.legendre.leggauss(npositive)
    mu = 0.5 * (nodes + 1.0)

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
        BDRF_Fourier_modes=arts.vdisort.fresnel_fourier_modes(
            refractive_index, 1
        ),
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


def test_fresnel_lambertian_depolarization():
    """Show how a depolarizing Lambertian component dilutes Fresnel Q."""
    refractive_index = 1.5
    lambertian_albedo = 0.2
    fresnel_fractions = (1.0, 0.5, 0.1, 0.02, 0.0)
    nquad = 24
    npositive = nquad // 2
    depth = 1.0e-8

    nodes, _ = np.polynomial.legendre.leggauss(npositive)
    mu = 0.5 * (nodes + 1.0)
    quadrature_angle = np.arccos(mu)
    phase = np.zeros((2, 1, 1, nquad, nquad, 4, 4))
    boundary_down = np.zeros((2, 1, npositive, 4))
    boundary_down[arts.vdisort.cosine_mode, 0, :, 0] = 1.0

    def reflected_stokes(modes):
        model = arts.cppvdisort(
            tau_arr=np.array([depth]),
            omega_arr=np.array([0.0]),
            NQuad=nquad,
            phase_matrix=phase,
            mu0=0.5,
            beam_stokes=np.zeros(4),
            phi0=0.0,
            b_neg=boundary_down,
            BDRF_Fourier_modes=modes,
        )
        return np.asarray(model.u(np.array([depth]), np.array([0.0])))[
            0, 0, :npositive
        ]

    pure_fresnel = reflected_stokes(
        arts.vdisort.fresnel_fourier_modes(refractive_index, 1)
    )
    pure_lambertian = reflected_stokes(
        arts.vdisort.lambertian_fourier_modes(lambertian_albedo, 1)
    )

    reflected = {}
    degree_of_linear_polarization = {}
    previous_peak = np.inf
    for fraction in fresnel_fractions:
        value = reflected_stokes(
            arts.vdisort.fresnel_lambertian_fourier_modes(
                fraction, lambertian_albedo, refractive_index, 1
            )
        )
        expected = fraction * pure_fresnel + (1.0 - fraction) * pure_lambertian
        np.testing.assert_allclose(value, expected, rtol=2.0e-8, atol=2.0e-12)
        np.testing.assert_allclose(value[:, 2:], 0.0, atol=2.0e-12)

        dolp = np.divide(
            np.abs(value[:, 1]),
            value[:, 0],
            out=np.zeros_like(value[:, 0]),
            where=value[:, 0] > 0.0,
        )
        reflected[fraction] = value
        degree_of_linear_polarization[fraction] = dolp
        assert np.max(dolp) <= previous_peak + 2.0e-12
        previous_peak = np.max(dolp)

    np.testing.assert_allclose(
        degree_of_linear_polarization[0.0], 0.0, atol=2.0e-12
    )
    assert np.max(degree_of_linear_polarization[1.0]) > 0.99
    assert np.max(degree_of_linear_polarization[0.02]) < 0.1

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        angle = np.linspace(0.0, 0.5 * np.pi - 1.0e-6, 500)
        dense_rv, dense_rh = _fresnel_amplitudes(
            np.cos(angle), refractive_index
        )
        fresnel_i = 0.5 * (dense_rv**2 + dense_rh**2)
        fresnel_q = 0.5 * (dense_rv**2 - dense_rh**2)

        _, (intensity_plot, polarization_plot) = plt.subplots(
            1, 2, figsize=(12, 4), constrained_layout=True
        )
        for fraction in fresnel_fractions:
            mixed_i = (
                fraction * fresnel_i
                + (1.0 - fraction) * lambertian_albedo
            )
            mixed_q = fraction * fresnel_q
            dense_dolp = np.divide(
                np.abs(mixed_q),
                mixed_i,
                out=np.zeros_like(mixed_i),
                where=mixed_i > 0.0,
            )
            label = f"{100.0 * fraction:g}% Fresnel"
            line = intensity_plot.plot(
                np.rad2deg(angle), mixed_i, label=label
            )[0]
            intensity_plot.plot(
                np.rad2deg(quadrature_angle),
                reflected[fraction][:, 0],
                "o",
                color=line.get_color(),
                fillstyle="none",
                markersize=4,
            )
            polarization_plot.plot(
                np.rad2deg(angle),
                100.0 * dense_dolp,
                color=line.get_color(),
                label=label,
            )
            polarization_plot.plot(
                np.rad2deg(quadrature_angle),
                100.0 * degree_of_linear_polarization[fraction],
                "o",
                color=line.get_color(),
                fillstyle="none",
                markersize=4,
            )

        brewster_angle = np.rad2deg(np.arctan(refractive_index))
        for axis in (intensity_plot, polarization_plot):
            axis.axvline(
                brewster_angle,
                color="black",
                linestyle="--",
                label="Brewster angle" if axis is intensity_plot else None,
            )
            axis.set_xlim(0.0, 90.0)
            axis.grid(True, alpha=0.25)
        intensity_plot.set(
            xlabel="Positive incidence angle [degree]",
            ylabel="Reflected Stokes I",
            title=f"Fresnel + Lambertian (A={lambertian_albedo:g})",
        )
        polarization_plot.set(
            xlabel="Positive incidence angle [degree]",
            ylabel="Degree of linear polarization [%]",
            title="Depolarization by the Lambertian component",
            ylim=(-2.0, 102.0),
        )
        intensity_plot.legend(ncols=2)
        polarization_plot.legend(ncols=2)
        plt.show()


def test_optically_thin_rayleigh_atmosphere():
    """Recover the analytical single-Rayleigh-scattering Stokes field."""
    depth = 1.0e-8
    nquad = 4
    nfourier = 1

    # A vertical beam makes the Rayleigh field axisymmetric.  Diffuse
    # rescattering is omitted from the transport phase matrix; its contribution
    # is second order in this optical depth and is outside the tested
    # single-scattering limit.
    phase = np.zeros((2, nfourier, 1, nquad, nquad, 4, 4))
    quadrature_beam_phase = np.zeros((2, nfourier, 1, nquad, 4, 4))
    quadrature_nodes, _ = np.polynomial.legendre.leggauss(nquad // 2)
    positive_mu = 0.5 * (quadrature_nodes + 1.0)
    signed_mu = np.concatenate((positive_mu, -positive_mu))
    quadrature_rayleigh = _rayleigh_phase_matrix(-signed_mu)
    quadrature_beam_phase[arts.vdisort.cosine_mode, 0, 0, :, :2, :2] = (
        quadrature_rayleigh[:, :2, :2]
    )
    quadrature_beam_phase[arts.vdisort.sine_mode, 0, 0, :, 2:, 2:] = (
        quadrature_rayleigh[:, 2:, 2:]
    )
    model = arts.cppvdisort(
        tau_arr=np.array([depth]),
        omega_arr=np.array([1.0]),
        NQuad=nquad,
        phase_matrix=phase,
        mu0=1.0,
        beam_stokes=np.array([1.0, 0.0, 0.0, 0.0]),
        phi0=0.0,
        beam_phase_matrix=quadrature_beam_phase,
    )

    # Positive mu observes the upward field at the top of the atmosphere.  For
    # a downward vertical beam, cos(scattering angle) = -mu.  The smallest mu
    # samples 90 degrees to within 0.006 degrees while remaining a valid ray.
    mu = np.geomspace(1.0e-4, 1.0, 80)
    rayleigh = _rayleigh_phase_matrix(-mu)
    user_phase = np.zeros((2, nfourier, 1, len(mu), nquad, 4, 4))
    user_beam_phase = np.zeros((2, nfourier, 1, len(mu), 4, 4))
    user_beam_phase[arts.vdisort.cosine_mode, 0, 0, :, :2, :2] = rayleigh[:, :2, :2]
    user_beam_phase[arts.vdisort.sine_mode, 0, 0, :, 2:, 2:] = rayleigh[:, 2:, 2:]

    stokes = np.asarray(
        model.u_user(
            tau=np.array([0.0]),
            phi=np.array([0.0]),
            mu=mu,
            phase_matrix=user_phase,
            beam_phase_matrix=user_beam_phase,
        )
    )[0, 0]

    # Along an upward ray at the top boundary, the attenuated direct-beam
    # source integrates to [1-exp(-depth*(1+1/mu))]/(1+mu).
    ray_integral = -np.expm1(-depth * (1.0 + 1.0 / mu)) / (1.0 + mu)
    expected = rayleigh[:, :, 0] * ray_integral[:, None] / (4.0 * np.pi)
    np.testing.assert_allclose(stokes, expected, rtol=2.0e-8, atol=2.0e-15)

    linear_polarization = np.hypot(stokes[:, 1], stokes[:, 2]) / stokes[:, 0]
    expected_polarization = (1.0 - mu**2) / (1.0 + mu**2)
    np.testing.assert_allclose(
        linear_polarization, expected_polarization, rtol=2.0e-8, atol=2.0e-12
    )
    assert linear_polarization[0] > 1.0 - 3.0e-8
    np.testing.assert_allclose(stokes[:, 2:], 0.0, atol=2.0e-15)

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        scattering_angle = np.rad2deg(np.arccos(-mu))
        _, (stokes_plot, polarization_plot) = plt.subplots(
            1, 2, figsize=(11, 4), constrained_layout=True
        )
        for component, label in enumerate("IQUV"):
            stokes_plot.plot(
                scattering_angle, stokes[:, component], label=f"VDISORT {label}"
            )
            stokes_plot.plot(
                scattering_angle,
                expected[:, component],
                "--",
                label=f"Analytical {label}",
            )
        stokes_plot.set(
            xlabel="Scattering angle [degree]",
            ylabel="Radiance",
            title="Optically thin Rayleigh Stokes field",
            xlim=(90.0, 180.0),
        )
        stokes_plot.grid(True, alpha=0.25)
        stokes_plot.legend(ncol=2)

        polarization_plot.plot(
            scattering_angle,
            linear_polarization,
            label="VDISORT",
        )
        polarization_plot.plot(
            scattering_angle,
            expected_polarization,
            "--",
            label="Analytical Rayleigh",
        )
        polarization_plot.set(
            xlabel="Scattering angle [degree]",
            ylabel="Degree of linear polarization",
            title="Rayleigh polarization",
            xlim=(90.0, 180.0),
            ylim=(-0.02, 1.02),
        )
        polarization_plot.grid(True, alpha=0.25)
        polarization_plot.legend()
        plt.show()


def test_nonprincipal_plane_rayleigh_polarization():
    """Reconstruct rotated Rayleigh Q/U outside the solar principal plane."""
    depth = 1.0e-8
    nquad = 4
    nfourier = 3
    mu0 = 0.6
    phi0 = 0.0

    phase = np.zeros((2, nfourier, 1, nquad, nquad, 4, 4))
    quadrature_nodes, _ = np.polynomial.legendre.leggauss(nquad // 2)
    positive_mu = 0.5 * (quadrature_nodes + 1.0)
    signed_mu = np.concatenate((positive_mu, -positive_mu))
    quadrature_beam_phase = _rayleigh_beam_fourier(signed_mu, mu0, phi0, nfourier)

    model = arts.cppvdisort(
        tau_arr=np.array([depth]),
        omega_arr=np.array([1.0]),
        NQuad=nquad,
        NFourier=nfourier,
        phase_matrix=phase,
        mu0=mu0,
        beam_stokes=np.array([1.0, 0.0, 0.0, 0.0]),
        phi0=phi0,
        beam_phase_matrix=quadrature_beam_phase,
    )

    user_mu = np.array([0.4])
    phi = np.linspace(0.0, 2.0 * np.pi, 180, endpoint=False)
    user_phase = np.zeros((2, nfourier, 1, 1, nquad, 4, 4))
    user_beam_phase = _rayleigh_beam_fourier(user_mu, mu0, phi0, nfourier)
    stokes = np.asarray(
        model.u_user(
            tau=np.array([0.0]),
            phi=phi,
            mu=user_mu,
            phase_matrix=user_phase,
            beam_phase_matrix=user_beam_phase,
        )
    )[0, :, 0]

    expected_phase = _rayleigh_stokes_column(user_mu[0], phi, mu0, phi0)
    ray_integral = -np.expm1(-depth * (1.0 / mu0 + 1.0 / user_mu[0])) / (
        1.0 + user_mu[0] / mu0
    )
    expected = expected_phase * ray_integral / (4.0 * np.pi)
    np.testing.assert_allclose(stokes, expected, rtol=3.0e-8, atol=3.0e-15)

    degree_linear = np.hypot(stokes[:, 1], stokes[:, 2]) / stokes[:, 0]
    cos_scattering = -user_mu[0] * mu0 + np.sqrt(1.0 - user_mu[0] ** 2) * np.sqrt(
        1.0 - mu0**2
    ) * np.cos(phi - phi0)
    expected_degree = (1.0 - cos_scattering**2) / (1.0 + cos_scattering**2)
    np.testing.assert_allclose(
        degree_linear, expected_degree, rtol=3.0e-8, atol=3.0e-12
    )
    np.testing.assert_allclose(stokes[:, 3], 0.0, atol=3.0e-15)

    principal = np.array([0, np.argmin(np.abs(phi - np.pi))])
    np.testing.assert_allclose(stokes[principal, 2], 0.0, atol=3.0e-15)
    nonprincipal = np.argmax(np.abs(stokes[:, 2]))
    assert np.abs(stokes[nonprincipal, 2]) > 0.2 * stokes[nonprincipal, 0]

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        scattering_angle = np.rad2deg(np.arccos(cos_scattering))
        figure, (stokes_plot, invariant_plot) = plt.subplots(
            1, 2, figsize=(11, 4), constrained_layout=True
        )
        stokes_plot.plot(np.rad2deg(phi), stokes[:, 1], label="VDISORT Q")
        stokes_plot.plot(np.rad2deg(phi), stokes[:, 2], label="VDISORT U")
        stokes_plot.plot(np.rad2deg(phi), expected[:, 1], "--", label="Analytical Q")
        stokes_plot.plot(np.rad2deg(phi), expected[:, 2], "--", label="Analytical U")
        stokes_plot.set(
            xlabel="Relative azimuth [degree]",
            ylabel="Polarized radiance",
            title="Non-principal-plane Rayleigh polarization",
            xlim=(0.0, 360.0),
        )
        stokes_plot.grid(True, alpha=0.25)
        stokes_plot.legend(ncol=2)

        invariant_plot.plot(scattering_angle, degree_linear, label="VDISORT")
        invariant_plot.plot(
            scattering_angle,
            expected_degree,
            "--",
            label="Analytical Rayleigh",
        )
        invariant_plot.set(
            xlabel="Scattering angle [degree]",
            ylabel="Degree of linear polarization",
            title="Rotation-invariant polarization",
            ylim=(-0.02, 1.02),
        )
        invariant_plot.grid(True, alpha=0.25)
        invariant_plot.legend()
        plt.show()


def test_optically_thin_sun_halo():
    """Transport a prescribed 22-degree halo through a thin atmosphere."""
    depth = 1.0e-8
    nquad = 4
    nfourier = 1

    # Normalize the prescribed phase function to unit spherical mean.  It is
    # axisymmetric about the vertical incident beam, so only Fourier mode zero
    # is present.
    normalization_nodes, normalization_weights = np.polynomial.legendre.leggauss(256)
    phase_normalization = 0.5 * np.dot(
        normalization_weights,
        _sun_halo_profile(np.arccos(normalization_nodes)),
    )

    phase = np.zeros((2, nfourier, 1, nquad, nquad, 4, 4))
    quadrature_beam_phase = np.zeros((2, nfourier, 1, nquad, 4, 4))
    quadrature_nodes, _ = np.polynomial.legendre.leggauss(nquad // 2)
    positive_mu = 0.5 * (quadrature_nodes + 1.0)
    signed_mu = np.concatenate((positive_mu, -positive_mu))
    quadrature_scattering_angle = np.arccos(-signed_mu)
    quadrature_beam_phase[arts.vdisort.cosine_mode, 0, 0, :, 0, 0] = (
        _sun_halo_profile(quadrature_scattering_angle) / phase_normalization
    )

    model = arts.cppvdisort(
        tau_arr=np.array([depth]),
        omega_arr=np.array([1.0]),
        NQuad=nquad,
        phase_matrix=phase,
        mu0=1.0,
        beam_stokes=np.array([1.0, 0.0, 0.0, 0.0]),
        phi0=0.0,
        beam_phase_matrix=quadrature_beam_phase,
    )

    scattering_angle = np.deg2rad(np.linspace(1.0, 45.0, 177))
    mu = -np.cos(scattering_angle)
    halo_phase = _sun_halo_profile(scattering_angle) / phase_normalization
    user_phase = np.zeros((2, nfourier, 1, len(mu), nquad, 4, 4))
    user_beam_phase = np.zeros((2, nfourier, 1, len(mu), 4, 4))
    user_beam_phase[arts.vdisort.cosine_mode, 0, 0, :, 0, 0] = halo_phase

    stokes = np.asarray(
        model.u_user(
            tau=np.array([depth]),
            phi=np.array([0.0]),
            mu=mu,
            phase_matrix=user_phase,
            beam_phase_matrix=user_beam_phase,
        )
    )[0, 0]

    abs_mu = np.abs(mu)
    slant_difference = depth * (1.0 / abs_mu - 1.0)
    ray_integral = np.exp(-depth) * -np.expm1(-slant_difference) / (1.0 - abs_mu)
    expected_intensity = halo_phase * ray_integral / (4.0 * np.pi)
    np.testing.assert_allclose(
        stokes[:, 0], expected_intensity, rtol=2.0e-8, atol=2.0e-15
    )
    np.testing.assert_allclose(stokes[:, 1:], 0.0, atol=2.0e-15)

    peak = np.argmax(stokes[:, 0])
    peak_angle = scattering_angle[peak]
    assert abs(peak_angle - np.deg2rad(22.0)) <= np.deg2rad(0.25)
    assert stokes[peak, 0] > 40.0 * stokes[0, 0]
    assert stokes[peak, 0] > 25.0 * stokes[-1, 0]

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        figure = plt.figure(figsize=(11, 4), constrained_layout=True)
        angular_plot = figure.add_subplot(1, 2, 1)
        sky_plot = figure.add_subplot(1, 2, 2, projection="polar")

        angular_plot.plot(np.rad2deg(scattering_angle), stokes[:, 0], label="VDISORT")
        angular_plot.plot(
            np.rad2deg(scattering_angle),
            expected_intensity,
            "--",
            label="Analytical single scattering",
        )
        angular_plot.axvline(22.0, color="black", linestyle=":")
        angular_plot.set(
            xlabel="Angular distance from the Sun [degree]",
            ylabel="Radiance",
            title="22-degree halo angular cut",
            xlim=(0.0, 45.0),
        )
        angular_plot.grid(True, alpha=0.25)
        angular_plot.legend()

        sky_azimuth = np.linspace(0.0, 2.0 * np.pi, 181)
        azimuth_grid, radius_grid = np.meshgrid(
            sky_azimuth, scattering_angle, indexing="xy"
        )
        halo_image = np.broadcast_to(stokes[:, 0, None], radius_grid.shape)
        image = sky_plot.pcolormesh(
            azimuth_grid,
            radius_grid,
            halo_image,
            shading="auto",
            cmap="inferno",
        )
        sky_plot.set(
            title="Axisymmetric sky view",
            rmax=np.deg2rad(45.0),
            rticks=np.deg2rad([10.0, 22.0, 30.0, 45.0]),
        )
        sky_plot.set_yticklabels(["10°", "22°", "30°", "45°"])
        figure.colorbar(image, ax=sky_plot, label="Radiance", pad=0.1)
        plt.show()


if __name__ == "__main__":
    test_combined_surface_models()
    test_phase_helpers()
    test_bdrf()
    test_absorbing_stokes_field()
    test_fresnel_brewster_angle()
    test_fresnel_lambertian_depolarization()
    test_optically_thin_rayleigh_atmosphere()
    test_nonprincipal_plane_rayleigh_polarization()
    test_optically_thin_sun_halo()
