#!/usr/bin/env python3

"""Smoke tests for the low-level polarized VDISORT Python interface."""

import os

import numpy as np

from pyarts3 import arts


def _fresnel_reflectances(mu, refractive_index):
    """VDISORT's Fresnel reflectances (vertical, horizontal) of unpolarized light from M00 +- M01."""
    matrices = np.array(
        [np.asarray(arts.vdisort.fresnel_reflection(x, refractive_index)) for x in np.atleast_1d(mu)]
    )
    return matrices[:, 0, 0] + matrices[:, 0, 1], matrices[:, 0, 0] - matrices[:, 0, 1]


def _rayleigh_beam_phase_matrix(signed_mu, mu0, nfourier):
    """ARTS's Rayleigh Fourier modes from the beam at -mu0 into signed_mu, in VDISORT's combined layout."""
    species = arts.ArrayOfScatteringSpecies(
        [arts.GasScatterer(arts.ConstantGasScattering(1e-30), arts.RayleighGasScattering(0.0))]
    )
    atm = arts.AtmPoint()
    atm.pressure = 8e4
    atm.temperature = 260.0
    optics = arts.vdisort.scattering_optics(species, atm, arts.AtmSpectralPoint(), 89e9, signed_mu, [-mu0], nfourier)
    return arts.vdisort.combine_beam_phase_matrices(
        np.asarray(optics.cosine)[:, None, :, 0], np.asarray(optics.sine)[:, None, :, 0]
    )


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
        lambertian_albedo, number_of_modes, mu, mu
    )
    fresnel = arts.vdisort.fresnel_fourier_modes(1.5, number_of_modes, mu, mu)
    generic = arts.vdisort.combine_fourier_modes(
        fresnel, fresnel_fraction, lambertian, 1.0 - fresnel_fraction
    )
    named = arts.vdisort.fresnel_lambertian_fourier_modes(
        fresnel_fraction, lambertian_albedo, 1.5, number_of_modes, mu, mu
    )

    # The kernels mix, and the specular part is the Fresnel one, weighted
    expected = fresnel_fraction * np.asarray(fresnel.brdf)
    expected += (1.0 - fresnel_fraction) * np.asarray(lambertian.brdf)
    for mixture in (generic, named):
        np.testing.assert_allclose(mixture.brdf, expected)
        np.testing.assert_allclose(
            mixture.specular, fresnel_fraction * np.asarray(fresnel.specular)
        )
    assert not lambertian.has_specular

    first_lambertian = np.asarray(lambertian.brdf)[arts.vdisort.cosine_mode, 0]
    np.testing.assert_allclose(
        first_lambertian[:, :, 0, 0], 2.0 * lambertian_albedo / np.pi
    )
    first_lambertian[:, :, 0, 0] = 0.0
    np.testing.assert_allclose(first_lambertian, 0.0)

    rough = arts.vdisort.cox_munk_lambertian_fourier_modes(
        0.35, lambertian_albedo, 5.0, 1.34, True, 2, mu, mu, azimuth_quadrature_points=24
    )
    assert rough.modes == 2
    assert np.isfinite(rough.brdf).all()

    scalar = arts.disort.cox_munk_lambertian_fourier_modes(
        0.35, lambertian_albedo, 5.0, 1.34, True, 2, mu, mu, azimuth_quadrature_points=24
    )
    assert scalar.modes == 2
    assert np.isfinite(scalar.brdf).all()

    absorbing = np.asarray(arts.vdisort.fresnel_reflection(0.5, 1.5 + 0.1j))
    assert np.isfinite(absorbing).all()
    assert abs(absorbing[2, 3]) > 1.0e-12


def test_delta_m_preprocessing():
    physical_tau = np.array([1.0, 3.0])
    physical_omega = np.array([0.5, 0.8])
    fraction = np.array([0.2, 0.25])
    phase = np.full((2, 1, 2, 2, 2, 4, 4), 3.0)
    removed = np.full_like(phase, 0.5)
    beam_phase = np.full((2, 1, 2, 2, 4, 4), 4.0)
    removed_beam = np.full_like(beam_phase, 1.0)

    transport = arts.vdisort.delta_m_preprocess(
        physical_tau,
        physical_omega,
        fraction,
        phase,
        removed,
        beam_phase,
        removed_beam,
    )
    np.testing.assert_allclose(transport.tau, [0.9, 2.5])
    np.testing.assert_allclose(transport.omega, [0.4444444444444444, 0.75])
    np.testing.assert_allclose(transport.source_coordinate_scale, [1.1111111111111112, 1.25])
    np.testing.assert_allclose(transport.source_coordinate_offset, [0.0, -0.125])
    np.testing.assert_allclose(np.asarray(transport.phase_matrix)[:, :, 0], 3.625)
    np.testing.assert_allclose(np.asarray(transport.phase_matrix)[:, :, 1], 3.8333333333333335)
    np.testing.assert_allclose(np.asarray(transport.beam_phase_matrix)[:, :, 0], 4.75)
    np.testing.assert_allclose(np.asarray(transport.beam_phase_matrix)[:, :, 1], 5.0)


def test_depolarizing_surface_catalogue():
    mu = np.array([0.2, 0.7])
    models = (
        arts.vdisort.hapke_fourier_modes(1.0, 0.06, 0.6, 2, mu, mu, azimuth_quadrature_points=32),
        arts.vdisort.rpv_fourier_modes(
            0.027, 0.647, -0.169, 0.1, 2, mu, mu, azimuth_quadrature_points=32
        ),
        arts.vdisort.ross_li_fourier_modes(
            0.091, 0.02, 0.01, np.deg2rad(1.5), 2, mu, mu, azimuth_quadrature_points=32
        ),
    )
    for model in models:
        assert model.modes == 2
        brdf = np.asarray(model.brdf)
        assert np.isfinite(brdf).all()
        # Fully depolarizing: only the I-to-I element of every Mueller block
        brdf[..., 0, 0] = 0.0
        np.testing.assert_allclose(brdf, 0.0, atol=1e-14)


def test_depolarizing_surface_catalogue_combinations():
    """Plot the new BRDFs alone and diluted by a polarized Fresnel component."""
    nquad = 16
    npositive = nquad // 2
    depth = 1.0e-8
    specular_fraction = 0.2
    refractive_index = 1.5

    mu = np.asarray(arts.math.pdleggauss(nquad)[0])  # VDISORT's upward streams
    viewing_angle = np.rad2deg(np.arccos(mu))
    phase = np.zeros((2, 1, 1, nquad, nquad, 4, 4))
    boundary_down = np.zeros((2, 1, npositive, 4))
    boundary_down[arts.vdisort.cosine_mode, 0, :, 0] = 1.0

    surfaces = {
        "Hapke": arts.vdisort.hapke_fourier_modes(
            1.0, 0.06, 0.6, 1, mu, mu, azimuth_quadrature_points=48
        ),
        "RPV": arts.vdisort.rpv_fourier_modes(
            0.027, 0.647, -0.169, 0.1, 1, mu, mu, azimuth_quadrature_points=48
        ),
        "Ross-Li": arts.vdisort.ross_li_fourier_modes(
            0.091, 0.02, 0.01, np.deg2rad(1.5), 1, mu, mu, azimuth_quadrature_points=48
        ),
    }
    fresnel = arts.vdisort.fresnel_fourier_modes(refractive_index, 1, mu, mu)

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
            surface=modes,
        )
        return np.asarray(model.u(np.array([depth]), np.array([0.0])))[0, 0, :npositive]

    fresnel_stokes = reflected_stokes(fresnel)
    pure = {}
    mixed = {}
    for name, surface in surfaces.items():
        pure[name] = reflected_stokes(surface)
        mixture = arts.vdisort.combine_fourier_modes(
            fresnel, specular_fraction, surface, 1.0 - specular_fraction
        )
        mixed[name] = reflected_stokes(mixture)
        np.testing.assert_allclose(pure[name][:, 1:], 0.0, atol=2.0e-12)
        np.testing.assert_allclose(
            mixed[name],
            specular_fraction * fresnel_stokes + (1.0 - specular_fraction) * pure[name],
            rtol=2.0e-8,
            atol=2.0e-12,
        )
        np.testing.assert_allclose(mixed[name][:, 2:], 0.0, atol=2.0e-12)

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        pure_intensity, mixed_intensity, mixed_q, polarization = axes.flat
        order = np.argsort(viewing_angle)
        angle = viewing_angle[order]

        for name in surfaces:
            line = pure_intensity.plot(angle, pure[name][order, 0], "o-", label=name)[0]
            color = line.get_color()
            mixed_intensity.plot(
                angle, mixed[name][order, 0], "o-", color=color, label=name
            )
            mixed_q.plot(angle, mixed[name][order, 1], "o-", color=color, label=name)
            degree = np.divide(
                np.abs(mixed[name][order, 1]),
                mixed[name][order, 0],
                out=np.zeros(npositive),
                where=mixed[name][order, 0] > 0.0,
            )
            polarization.plot(angle, 100.0 * degree, "o-", color=color, label=name)

        pure_intensity.set(
            title="Fully depolarizing BRDF catalogue",
            ylabel="Reflected Stokes I",
        )
        mixed_intensity.set(
            title=f"{100 * specular_fraction:g}% Fresnel + BRDF",
            ylabel="Reflected Stokes I",
        )
        mixed_q.set(
            title="Polarized contribution after mixing",
            xlabel="Positive viewing angle [degree]",
            ylabel="Reflected Stokes Q",
        )
        polarization.set(
            title="Depolarization of Fresnel reflection",
            xlabel="Positive viewing angle [degree]",
            ylabel="Degree of linear polarization [%]",
        )
        for axis in axes.flat:
            axis.set_xlim(0.0, 90.0)
            axis.grid(True, alpha=0.25)
            axis.legend()
        figure.suptitle("VDISORT surface BRDFs under isotropic illumination")
        plt.show()


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


def test_surface_reflection():
    """A surface reflection has the shapes it is constructed with; its data is set in place."""
    surface = arts.vdisort.SurfaceReflection(2, 3, 4, 0.5)
    assert (surface.modes, surface.directions, surface.streams) == (2, 3, 4)
    assert surface.beam_mu == 0.5 and not surface.has_specular
    np.testing.assert_array_equal(surface.brdf, np.zeros((2, 2, 3, 4, 4, 4)))

    brdf = np.ones((2, 2, 3, 4, 4, 4))
    beam_brdf = np.full((2, 2, 3, 4, 4), 2.0)
    surface.set(brdf, beam_brdf, 0.6)
    np.testing.assert_array_equal(surface.brdf, brdf)
    np.testing.assert_array_equal(surface.beam_brdf, beam_brdf)
    assert surface.beam_mu == 0.6

    # The reflection into other directions for the same modes, streams and beam
    user = surface.like(5)
    assert (user.modes, user.directions, user.streams, user.beam_mu) == (2, 5, 4, 0.6)
    np.testing.assert_array_equal(user.brdf, np.zeros((2, 2, 5, 4, 4, 4)))

    from_data = arts.vdisort.SurfaceReflection(brdf, beam_brdf, 0.6)
    np.testing.assert_array_equal(from_data.brdf, brdf)

    # Other shapes are refused, and a specular part cannot reflect a beam
    for refused in (
        lambda: surface.set(np.ones((2, 2, 3, 3, 4, 4)), beam_brdf, 0.6),
        lambda: surface.set(brdf, beam_brdf, 1.5),
        lambda: arts.vdisort.SurfaceReflection(1, 2, 2, 0.5, True),
    ):
        try:
            refused()
        except RuntimeError:
            pass
        else:
            raise AssertionError("A surface reflection accepted invalid data")
    np.testing.assert_array_equal(surface.brdf, brdf)


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
    mu = np.asarray(arts.math.pdleggauss(nquad)[0])

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
        surface=arts.vdisort.fresnel_fourier_modes(refractive_index, 1, mu, mu),
    )

    stokes = np.asarray(model.u(np.array([depth]), np.array([0.0])))[0, 0, :npositive]
    # VDISORT's own Fresnel matrix, attenuated by the transparent layer
    fresnel_vertical, fresnel_horizontal = _fresnel_reflectances(mu, refractive_index)
    attenuation = np.exp(-depth / mu)
    expected_vertical = fresnel_vertical * attenuation
    expected_horizontal = fresnel_horizontal * attenuation
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
        dense_vertical, dense_horizontal = _fresnel_reflectances(np.cos(angle), refractive_index)
        quadrature_angle = np.arccos(mu)

        _, (reflectance_plot, error_plot) = plt.subplots(
            1, 2, figsize=(11, 4), constrained_layout=True
        )
        reflectance_plot.plot(np.rad2deg(angle), dense_vertical, label="Fresnel p")
        reflectance_plot.plot(np.rad2deg(angle), dense_horizontal, label="Fresnel s")
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
            title="VDISORT minus its Fresnel matrix",
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

    mu = np.asarray(arts.math.pdleggauss(nquad)[0])  # VDISORT's upward streams
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
            surface=modes,
        )
        return np.asarray(model.u(np.array([depth]), np.array([0.0])))[0, 0, :npositive]

    pure_fresnel = reflected_stokes(
        arts.vdisort.fresnel_fourier_modes(refractive_index, 1, mu, mu)
    )
    pure_lambertian = reflected_stokes(
        arts.vdisort.lambertian_fourier_modes(lambertian_albedo, 1, mu, mu)
    )

    reflected = {}
    degree_of_linear_polarization = {}
    previous_peak = np.inf
    for fraction in fresnel_fractions:
        value = reflected_stokes(
            arts.vdisort.fresnel_lambertian_fourier_modes(
                fraction, lambertian_albedo, refractive_index, 1, mu, mu
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

    np.testing.assert_allclose(degree_of_linear_polarization[0.0], 0.0, atol=2.0e-12)
    assert np.max(degree_of_linear_polarization[1.0]) > 0.99
    assert np.max(degree_of_linear_polarization[0.02]) < 0.1

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        angle = np.linspace(0.0, 0.5 * np.pi - 1.0e-6, 500)
        dense_vertical, dense_horizontal = _fresnel_reflectances(np.cos(angle), refractive_index)
        fresnel_i = 0.5 * (dense_vertical + dense_horizontal)
        fresnel_q = 0.5 * (dense_vertical - dense_horizontal)

        _, (intensity_plot, polarization_plot) = plt.subplots(
            1, 2, figsize=(12, 4), constrained_layout=True
        )
        for fraction in fresnel_fractions:
            mixed_i = fraction * fresnel_i + (1.0 - fraction) * lambertian_albedo
            mixed_q = fraction * fresnel_q
            dense_dolp = np.divide(
                np.abs(mixed_q),
                mixed_i,
                out=np.zeros_like(mixed_i),
                where=mixed_i > 0.0,
            )
            label = f"{100.0 * fraction:g}% Fresnel"
            line = intensity_plot.plot(np.rad2deg(angle), mixed_i, label=label)[0]
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
    """Recover the single-Rayleigh-scattering Stokes field of a vertical beam."""
    depth = 1.0e-8
    nquad = 160
    nfourier = 1

    # A vertical beam makes the Rayleigh field axisymmetric.  Diffuse
    # rescattering is omitted from the transport phase matrix; its contribution
    # is second order in this optical depth and is outside the tested
    # single-scattering limit.
    phase = np.zeros((2, nfourier, 1, nquad, nquad, 4, 4))
    positive_mu, _ = arts.math.pdleggauss(nquad)
    signed_mu = np.concatenate((positive_mu, -positive_mu))
    model = arts.cppvdisort(
        tau_arr=np.array([depth]),
        omega_arr=np.array([1.0]),
        NQuad=nquad,
        phase_matrix=phase,
        mu0=1.0,
        beam_stokes=np.array([1.0, 0.0, 0.0, 0.0]),
        phi0=0.0,
        beam_phase_matrix=_rayleigh_beam_phase_matrix(signed_mu, 1.0, nfourier),
    )

    # The upward streams at the top of the atmosphere.  For a downward
    # vertical beam, cos(scattering angle) = -mu.  The smallest of the 80
    # streams samples 90 degrees to within 0.02 degrees.
    mu = np.asarray(model.mu)[: nquad // 2]
    stokes = np.asarray(model.u(tau=np.array([0.0]), phi=np.array([0.0])))[0, 0, : nquad // 2]

    # The single-scattering I and Q at every tenth stream and the last one
    streams = [0, 10, 20, 30, 40, 50, 60, 70, 79]
    expected = np.array([
        [2.675248706155e-06, -2.675248439868e-06],
        [1.378775389816e-08, -1.373598663198e-08],
        [3.938644033591e-09, -3.753410154038e-09],
        [2.062135924119e-09, -1.681539616990e-09],
        [1.475055794402e-09, -8.665797677233e-10],
        [1.270897602141e-09, -4.362930262823e-10],
        [1.207542718395e-09, -1.825658856973e-10],
        [1.194454033628e-09, -4.348933782705e-11],
        [1.193662090961e-09, -2.663221967161e-13],
    ])  # fmt: skip
    np.testing.assert_allclose(stokes[streams, :2], expected, rtol=2.0e-8, atol=2.0e-15)
    np.testing.assert_allclose(stokes[:, 2:], 0.0, atol=2.0e-15)

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        scattering_angle = np.rad2deg(np.arccos(-mu))
        _, (stokes_plot, polarization_plot) = plt.subplots(
            1, 2, figsize=(11, 4), constrained_layout=True
        )
        for component, label in enumerate("IQ"):
            line = stokes_plot.plot(
                scattering_angle, stokes[:, component], label=f"VDISORT {label}"
            )[0]
            stokes_plot.plot(
                scattering_angle[streams],
                expected[:, component],
                "o",
                color=line.get_color(),
                fillstyle="none",
                label=f"Expected {label}",
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
            np.hypot(stokes[:, 1], stokes[:, 2]) / stokes[:, 0],
            label="VDISORT",
        )
        polarization_plot.plot(
            scattering_angle[streams],
            np.abs(expected[:, 1]) / expected[:, 0],
            "o",
            fillstyle="none",
            label="Expected",
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
    """Rotated Rayleigh Q/U outside the solar principal plane."""
    depth = 1.0e-8
    nquad = 4
    nfourier = 3
    mu0 = 0.6
    phi0 = 0.0

    phase = np.zeros((2, nfourier, 1, nquad, nquad, 4, 4))
    positive_mu, _ = arts.math.pdleggauss(nquad)
    signed_mu = np.concatenate((positive_mu, -positive_mu))

    model = arts.cppvdisort(
        tau_arr=np.array([depth]),
        omega_arr=np.array([1.0]),
        NQuad=nquad,
        NFourier=nfourier,
        phase_matrix=phase,
        mu0=mu0,
        beam_stokes=np.array([1.0, 0.0, 0.0, 0.0]),
        phi0=phi0,
        beam_phase_matrix=_rayleigh_beam_phase_matrix(signed_mu, mu0, nfourier),
    )

    # The upward stream closest to mu = 0.4 (0.2113), at the top, every 45 degrees
    upward = np.asarray(model.mu)[: nquad // 2]
    stream = int(np.argmin(np.abs(upward - 0.4)))
    phi = np.deg2rad(np.arange(0.0, 360.0, 45.0))
    stokes = np.asarray(model.u(tau=np.array([0.0]), phi=phi))[0, :, stream]

    # The single-scattering I, Q and U
    expected = np.array([
        [4.036411810380e-09, -1.612057586305e-09, 0.0],
        [3.337042265115e-09, -5.039169246308e-10, 2.255828698644e-09],
        [2.869639782742e-09, 8.361907999352e-10, 2.650033767602e-09],
        [4.129027357264e-09, 2.880681675188e-10, 1.491884996245e-09],
        [5.156447868895e-09, -4.920215277896e-10, 0.0],
        [4.129027357264e-09, 2.880681675188e-10, -1.491884996245e-09],
        [2.869639782742e-09, 8.361907999352e-10, -2.650033767602e-09],
        [3.337042265115e-09, -5.039169246308e-10, -2.255828698644e-09],
    ])  # fmt: skip
    np.testing.assert_allclose(stokes[:, :3], expected, rtol=3.0e-8, atol=3.0e-15)
    np.testing.assert_allclose(stokes[:, 3], 0.0, atol=3.0e-15)

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        dense_phi = np.linspace(0.0, 2.0 * np.pi, 180, endpoint=False)
        dense = np.asarray(model.u(tau=np.array([0.0]), phi=dense_phi))[0, :, stream]
        figure, (stokes_plot, invariant_plot) = plt.subplots(
            1, 2, figsize=(11, 4), constrained_layout=True
        )
        for component, label in ((1, "Q"), (2, "U")):
            line = stokes_plot.plot(np.rad2deg(dense_phi), dense[:, component], label=f"VDISORT {label}")[0]
            stokes_plot.plot(
                np.rad2deg(phi), expected[:, component], "o", color=line.get_color(), fillstyle="none",
                label=f"Expected {label}",
            )
        stokes_plot.set(
            xlabel="Relative azimuth [degree]",
            ylabel="Polarized radiance",
            title="Non-principal-plane Rayleigh polarization",
            xlim=(0.0, 360.0),
        )
        stokes_plot.grid(True, alpha=0.25)
        stokes_plot.legend(ncol=2)

        invariant_plot.plot(
            np.rad2deg(dense_phi), np.hypot(dense[:, 1], dense[:, 2]) / dense[:, 0], label="VDISORT"
        )
        invariant_plot.plot(
            np.rad2deg(phi),
            np.hypot(expected[:, 1], expected[:, 2]) / expected[:, 0],
            "o",
            fillstyle="none",
            label="Expected",
        )
        invariant_plot.set(
            xlabel="Relative azimuth [degree]",
            ylabel="Degree of linear polarization",
            title="Degree of linear polarization",
            ylim=(-0.02, 1.02),
        )
        invariant_plot.grid(True, alpha=0.25)
        invariant_plot.legend()
        plt.show()


def test_optically_thin_sun_halo():
    """Transport a prescribed 22-degree halo through a thin atmosphere."""
    # Without diffuse scattering the single-scattering radiance is exact at
    # any depth.  The stream solution subtracts exp(-tau/mu) from
    # exp(-tau/mu0) with mu0 = 1, so depth 1e-3 rather than 1e-8 keeps its
    # cancellation (eps / (depth (1/mu - 1))) below the tolerance.
    depth = 1.0e-3
    nquad = 256
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
    positive_mu, _ = arts.math.pdleggauss(nquad)
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

    # The downward streams at the bottom with scattering angles from 1 to 45
    # degrees, in ascending scattering angle
    downward = np.asarray(model.mu)[nquad // 2 :]
    streams = np.flatnonzero(
        (-downward <= np.cos(np.deg2rad(1.0))) & (-downward >= np.cos(np.deg2rad(45.0)))
    )
    streams = streams[np.argsort(np.arccos(-downward[streams]))]
    mu = downward[streams]
    scattering_angle = np.arccos(-mu)
    stokes = np.asarray(model.u(tau=np.array([depth]), phi=np.array([0.0])))[0, 0, nquad // 2 + streams]

    # The single-scattering I of every third of the 46 streams, 1.7 to 45 degrees
    expected_intensity = np.array([
        4.371866781738e-05, 4.384626307824e-05, 4.409303287014e-05, 4.446166199835e-05,
        4.495620072908e-05, 4.583830361760e-05, 5.221255335508e-04, 3.650890152195e-03,
        1.690112990286e-04, 4.958999024698e-05, 5.099909306138e-05, 5.263398660074e-05,
        5.449749851022e-05, 5.661651543834e-05, 5.902317522057e-05, 6.175599818999e-05,
    ])  # fmt: skip
    assert len(mu) == 46
    np.testing.assert_allclose(
        stokes[::3, 0], expected_intensity, rtol=2.0e-8, atol=2.0e-15
    )
    np.testing.assert_allclose(stokes[:, 1:], 0.0, atol=2.0e-15)

    # The brightest stream is the one closest to 22 degrees
    peak = np.argmax(stokes[:, 0])
    peak_angle = scattering_angle[peak]
    spacing = np.max(np.abs(np.diff(scattering_angle[max(peak - 1, 0) : peak + 2])))
    assert abs(peak_angle - np.deg2rad(22.0)) <= 0.5 * spacing + 1e-12
    assert stokes[peak, 0] > 40.0 * stokes[0, 0]
    assert stokes[peak, 0] > 25.0 * stokes[-1, 0]

    if "ARTS_HEADLESS" not in os.environ:
        import matplotlib.pyplot as plt

        figure = plt.figure(figsize=(11, 4), constrained_layout=True)
        angular_plot = figure.add_subplot(1, 2, 1)
        sky_plot = figure.add_subplot(1, 2, 2, projection="polar")

        angular_plot.plot(np.rad2deg(scattering_angle), stokes[:, 0], label="VDISORT")
        angular_plot.plot(
            np.rad2deg(scattering_angle[::3]),
            expected_intensity,
            "o",
            fillstyle="none",
            label="Expected single scattering",
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
    test_delta_m_preprocessing()
    test_depolarizing_surface_catalogue()
    test_depolarizing_surface_catalogue_combinations()
    test_phase_helpers()
    test_surface_reflection()
    test_absorbing_stokes_field()
    test_fresnel_brewster_angle()
    test_fresnel_lambertian_depolarization()
    test_optically_thin_rayleigh_atmosphere()
    test_nonprincipal_plane_rayleigh_polarization()
    test_optically_thin_sun_halo()
