#!/usr/bin/env python3

"""Smoke tests for the workspace-level scalar DISORT BRDF methods.

The numerical BRDF and Fourier reconstruction tests live in the C++ DISORT
test suite.  This test verifies that ARTS users can configure every model
through a workspace: the settings hold the surface's reflection between the
streams, and of the sun's beam at the sun's zenith angle, which must therefore
be set first.
"""

import numpy as np

from pyarts3 import arts
from pyarts3.workspace import Workspace


NQUAD = 8
NFOURIER = 32
ZENITH = 60.0


def configure(method, **kwargs):
    ws = Workspace()
    settings = arts.DisortSettings()
    settings.freq_grid = np.array([1.0, 2.0])
    settings.quadrature_dimension = NQUAD
    settings.fourier_mode_dimension = NFOURIER
    # The surface reflects the sun, so the sun is set first
    settings.solar_zenith_angle = np.full(2, ZENITH)
    ws.disort_settings = settings
    getattr(ws, method)(azimuth_quadrature_points=128, **kwargs)

    settings = ws.disort_settings
    n = NQUAD // 2
    assert np.asarray(settings.surface_brdf).shape == (2, NFOURIER, n, n)
    assert np.asarray(settings.surface_beam_brdf).shape == (2, NFOURIER, n)
    assert np.isfinite(settings.surface_brdf).all()
    np.testing.assert_allclose(settings.surface_beam_mu, np.cos(np.deg2rad(ZENITH)))
    assert np.asarray(settings.surface_beam_brdf).any()
    return ws


def test_models():
    configure("disort_settingsSurfaceHapke")
    configure(
        "disort_settingsSurfaceHapkeSpectral",
        opposition_amplitude=[0.9, 1.1],
        opposition_width=[0.05, 0.07],
        single_scattering_albedo=[0.5, 0.7],
    )
    configure("disort_settingsSurfaceCoxMunk", shadowing=0)
    configure("disort_settingsSurfaceCoxMunk", shadowing=1)
    configure("disort_settingsSurfaceCoxMunkSpectral", refractive_index=[1.33, 1.35])
    configure("disort_settingsSurfaceRPV")
    configure(
        "disort_settingsSurfaceRPVSpectral",
        rho0=[0.02, 0.04],
        kappa=[0.6, 0.7],
        asymmetry=[-0.2, -0.1],
        hotspot=[0.08, 0.12],
    )
    configure("disort_settingsSurfaceRossLi")
    configure(
        "disort_settingsSurfaceRossLiSpectral",
        isotropic=[0.08, 0.1],
        volumetric=[0.02, 0.02],
        geometric=[0.005, 0.015],
    )


def test_spectral_length_check():
    try:
        configure(
            "disort_settingsSurfaceHapkeSpectral",
            opposition_amplitude=[1.0, 1.0],
            opposition_width=[0.06, 0.06],
            single_scattering_albedo=[0.5],
        )
    except RuntimeError:
        pass
    else:
        raise AssertionError("A mismatched spectral BRDF parameter was accepted")


def test_sun_required():
    ws = Workspace()
    settings = arts.DisortSettings()
    settings.freq_grid = np.array([1.0, 2.0])
    ws.disort_settings = settings
    try:
        ws.disort_settingsSurfaceHapke()
    except RuntimeError:
        pass
    else:
        raise AssertionError("A surface was set without the sun")


if __name__ == "__main__":
    test_models()
    test_spectral_length_check()
    test_sun_required()
