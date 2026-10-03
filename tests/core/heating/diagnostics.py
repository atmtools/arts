"""Physics and array-contract checks independent of the ARTS2 adapters.

This file independently checks the supported formulas with a quadratic
flux on nonuniform pressure levels (both orderings and both boundaries),
isotropic hemispheres, frequency integration, analytic beam heating and
thermal cooling in both solvers, workspace DFDT sampling locations, spectral
extinction weighting, and invalid inputs. This separates preserved historical
behavior from the public scientific interface.
"""

import numpy as np
from pyarts3 import arts, Workspace
from pyarts3.recipe import heating_rates as hr


def test_pressure_derivative():
    pressure = np.array([100000.0, 75000.0, 30000.0, 10000.0])
    flux = 2e-8 * pressure**2 - 0.003 * pressure + 200.0
    cp = np.array([1000.0, 1001.0, 1002.0, 1003.0])
    gravity = np.array([9.8, 9.7, 9.6, 9.5])
    exact = (4e-8 * pressure - 0.003) * gravity / cp
    np.testing.assert_allclose(
        hr.from_flux(flux, pressure, cp, gravity), exact, atol=1e-18
    )
    np.testing.assert_allclose(
        hr.from_flux(flux[::-1], pressure[::-1], cp[::-1], gravity[::-1]),
        exact[::-1],
        atol=1e-18,
    )
    np.testing.assert_allclose(
        hr.from_flux(
            np.stack([flux, 2 * flux], axis=1),
            pressure,
            cp[:, None],
            gravity[:, None],
            axis=0,
        ),
        np.stack([exact, 2 * exact], axis=1),
        atol=1e-18,
    )
    np.testing.assert_allclose(
        hr.from_flux(np.ones(4), pressure, 1000.0, 9.8), 0.0, atol=1e-18
    )


def test_radiance_integration():
    nodes, weights = np.polynomial.legendre.leggauss(3)
    mu = np.r_[(nodes + 1) / 2, -(nodes + 1) / 2]
    weights = np.tile(weights / 2, 2)
    zenith = np.rad2deg(np.arccos(mu))
    up, down = hr.flux_from_radiance(np.where(mu < 0, 3.0, 2.0), zenith, weights)
    np.testing.assert_allclose([up, down], np.pi * np.array([3.0, 2.0]), rtol=1e-14)
    # Nonuniform frequency grid: a linear spectrum integrates exactly.
    frequency = np.array([1e11, 2e11, 5e11])
    spectrum = frequency[:, None] * np.ones((3, 6)) * 1e-12
    integral = 0.5e-12 * (frequency[-1] ** 2 - frequency[0] ** 2)
    np.testing.assert_allclose(hr.integrate_spectral(spectrum, frequency), integral)
    up, down = hr.flux_from_radiance(spectrum.T, zenith, weights, axis=0)
    np.testing.assert_allclose(up, np.pi * frequency * 1e-12)
    np.testing.assert_allclose(up, down)
    np.testing.assert_allclose(
        hr.integrate_spectral([2.0, 4.0], [1.0, 3.0], weights=[5.0, 7.0]), 38.0
    )
    np.testing.assert_allclose(hr.integrate_spectral([2.0], [1.0], weights=[5.0]), 10.0)


def test_dfdt():
    # A downward beam in a homogeneous absorbing layer: Fnet=-mu0*I0*exp(-tau/mu0),
    # DFDT=I0*exp(-tau/mu0). This fixes both sign and optical-depth normalization.
    tau = np.array([0.0, 0.2, 0.8])
    mu0, beam = 0.7, 2.0
    scalar = arts.cppdisort(
        [1.0], [0.0], 6, [[1.0]], mu0, beam, 0.0, NLeg=1, NFourier=1
    )
    vector = arts.cppvdisort(
        [1.0],
        [0.0],
        6,
        np.zeros((2, 1, 1, 6, 6, 4, 4)),
        mu0,
        [beam, 0.0, 0.0, 0.0],
        0.0,
        NFourier=1,
    )
    expected = beam * np.exp(-tau / mu0)
    for solver in (scalar, vector):
        flux = np.asarray(solver.flux(tau))
        np.testing.assert_allclose(flux[3], expected, rtol=1e-13)
        np.testing.assert_allclose(
            hr.from_optical_depth_derivative(flux[3], 0.002, 1.2, 1000.0),
            expected * 0.002 / 1200.0,
            rtol=1e-13,
        )

    # Emitting gas with cold boundaries must cool. Its analytic mean intensity
    # is obtained from the same double-Gauss angular quadrature, independently
    # of the transfer solver and heating recipe.
    nodes, weights = np.polynomial.legendre.leggauss(3)
    mu, weights = (nodes + 1) / 2, weights / 2
    b = 2.0
    source = np.zeros((1, 1, 4))
    source[0, 0, 0] = b
    scalar = arts.cppdisort(
        [1.0], [0.0], 6, [[1.0]], 0.0, 0.0, 0.0, NLeg=1, NFourier=1, s_poly_coeffs=[[b]]
    )
    vector = arts.cppvdisort(
        [1.0],
        [0.0],
        6,
        np.zeros((2, 1, 1, 6, 6, 4, 4)),
        0.0,
        [0.0, 0.0, 0.0, 0.0],
        0.0,
        NFourier=1,
        s_poly_coeffs=source,
    )
    expected = (
        -2
        * np.pi
        * b
        * ((np.exp(-tau[:, None] / mu) + np.exp(-(1 - tau[:, None]) / mu)) @ weights)
    )
    for solver in (scalar, vector):
        np.testing.assert_allclose(
            np.asarray(solver.flux(tau))[3], expected, rtol=1e-13
        )


def test_workspace_dfdt():
    points = []
    for altitude in (1000.0, 500.0, 0.0):
        point = arts.PropagationPathPoint()
        point.pos = [altitude, 0.0, 0.0]
        points.append(point)
    ws = Workspace()
    ws.disort_settingsInit(
        freq_grid=[1e11, 2e11],
        ray_path=points,
        disort_quadrature_dimension=6,
        disort_legendre_polynomial_dimension=1,
        disort_fourier_mode_dimension=1,
    )
    s = ws.disort_settings
    s.optical_thicknesses = [[0.2, 0.6], [0.4, 1.2]]
    s.legendre_coefficients = np.ones((2, 2, 1))
    s.solar_source = [2.0, 3.0]
    s.solar_zenith_angle = [60.0, 60.0]
    ws.disort_spectral_flux_fieldCalc()
    flux = ws.disort_spectral_flux_field
    # Gridded output is at optical_thicknesses, not midpoints.
    expected = np.array([2.0, 3.0])[:, None] * np.exp(
        -np.asarray(s.optical_thicknesses) / 0.5
    )
    np.testing.assert_allclose(flux.dfdt, expected, rtol=1e-13)
    extinction = np.array([[0.2, 0.4], [0.4, 0.8]]) / 500.0
    density = [0.8, 1.2]
    expected_heat = np.sum(
        expected * extinction * np.array([2.0, 5.0])[:, None], axis=0
    ) / (np.array(density) * 1000.0)
    np.testing.assert_allclose(
        hr.from_disort(flux, extinction, density, 1000.0, weights=[2.0, 5.0]),
        expected_heat,
        rtol=1e-13,
    )


def test_invalid_inputs():
    for pressure in ([1.0, 1.0, 2.0], [1.0, 3.0, 2.0], [1.0, 2.0], [0.0, 1.0, 2.0]):
        with np.testing.assert_raises(ValueError):
            hr.from_flux(np.ones(len(pressure)), pressure, 1000.0, 9.8)
    for capacity in (0.0, -1.0, np.nan):
        with np.testing.assert_raises(ValueError):
            hr.from_flux_divergence(1.0, capacity, 9.8)
    with np.testing.assert_raises(ValueError):
        hr.from_optical_depth_derivative(1.0, -1.0, 1.0, 1000.0)
    with np.testing.assert_raises(ValueError):
        hr.from_optical_depth_derivative(1.0, 1.0, 0.0, 1000.0)
    for frequency in ([1.0], [2.0, 1.0], [1.0, np.nan]):
        with np.testing.assert_raises(ValueError):
            hr.integrate_spectral(np.ones(len(frequency)), frequency)
    with np.testing.assert_raises(ValueError):
        hr.integrate_spectral([1.0, 2.0], [1.0, 2.0], weights=[1.0])
    with np.testing.assert_raises(ValueError):
        hr.flux_from_radiance([1.0, 2.0], [0.0, 190.0], [1.0, 1.0])


test_pressure_derivative()
test_radiance_integration()
test_dfdt()
test_workspace_dfdt()
test_invalid_inputs()
