#!/usr/bin/env python3

"""Spectral delta-M-plus configuration through the ARTS workspace."""

import numpy as np

from pyarts3 import arts
from pyarts3.workspace import Workspace


def test_spectral_delta_m_plus():
    coefficients = np.array(
        [
            [[1.0, 0.9, 0.8, 0.72]],
            [[1.0, 0.85, 0.7, 0.60]],
        ]
    )

    settings = arts.DisortSettings()
    settings.freq_grid = np.array([1.0, 2.0])
    settings.alt_grid = np.array([1.0, 0.0])
    settings.legendre_polynomial_dimension = 2
    settings.legendre_coefficients = coefficients
    settings.fractional_scattering = np.zeros((2, 1))

    ws = Workspace()
    ws.disort_settings = settings
    ws.disort_settingsDeltaMPlus()

    actual_fraction = np.array(ws.disort_settings.fractional_scattering)
    actual_moments = np.array(ws.disort_settings.delta_m_peak_moments)
    for iv in range(2):
        expected_fraction, expected_moments = arts.disort.delta_m_plus(
            coefficients[iv], 2
        )
        np.testing.assert_allclose(actual_fraction[iv], expected_fraction)
        np.testing.assert_allclose(actual_moments[iv], expected_moments)


def test_standard_agenda_option():
    ws = Workspace()
    ws.disort_settings_agendaSetup(
        scattering_setting="ScatteringSpeciesDeltaMPlus"
    )


if __name__ == "__main__":
    test_spectral_delta_m_plus()
    test_standard_agenda_option()
