# -*- coding: utf-8 -*-

"""Functions directly related to spectroscopy.
"""
import numpy as np
from scipy import interpolate

import pyarts.constants as constants


__all__ = [
    'linewidth',
    'doppler_broadening',
    'boltzmann_level',
    'stimulated_emission',
]


def linewidth(f, a):
    """Calculate the full-width at half maximum (FWHM) of an absorption line.

    Parameters:
        f (ndarray): Frequency grid.
        a (ndarray): Line properties
            (e.g. absorption coefficients or cross-sections).

    Returns:
        float: Linewidth.

    Examples:
        >>> f = np.linspace(0, np.pi, 100)
        >>> a = np.sin(f)**2
        >>> linewidth(f, a)
        1.571048056449009
    """
    s = interpolate.UnivariateSpline(f, a - np.max(a)/2, s=0)
    return float(np.diff(s.roots()))


def doppler_broadening(t, f0, m):
    """Calculate the doppler broadening half-width half-maximum

    .. math::
        \\gamma_D(T) = \\sqrt{ \\frac{2\\log(2) k_B T}{mc^2} } f_0

    Parameters:
        t (float or ndarray): Temperature [Kelvin]

        f0 (float or like temperature): Central frequency [Hertz/invcm]

        m (float or like temperature): Mass [kilogram]

    Returns
        hwhm (like temperature): Half-width half-maximum [Hertz/invcm]
    """

    return np.sqrt(2 * constants.boltzmann * t * np.log(2) /
                   (m * constants.speed_of_light**2)) * f0


def boltzmann_level(elow, t, t0):
    """Computes the Boltzmann level function

    .. math::
        K_1 = \\exp\\left(\\frac{E_l \\left[T-T_0\\right]}{k_B T T_0}\\right),

    where :math:`k_B` is the Boltzmann constant.

    All ndarrays must be of same size, any of the inputs can be ndarray

    Parameters:
        elow (float or ndarray): Lower state energy level [J]

        t (float or ndarray): Temperature [Kelvin]

        t0 (float or ndarray): Line temperature [Kelvin]

    Returns
        K1 (like input): How much Boltzmann statistics feeds the transition

    .. math::
        S(T) = S(T_0)K_1K_2 \\frac{Q(T_0)}{Q(T)}
    """
    return np.exp(elow * (t - t0) / (constants.boltzmann * t * t0))


def stimulated_emission(f0, t, t0):
    """Computes the stimulated emission function

    .. math::
        K_2 = \\frac{1 - \\exp\\left( - \\frac{h f_0}{k_B T}\\right)}
        {1 - \\exp\\left( - \\frac{h f_0}{k_B T_0}\\right)},

    with Planck constant :math:`h` and Boltzmann constant :math:`k_B`.

    Parameters:
        f0 (float or ndarray): Line frequency [Hz]

        t (float or ndarray): Temperature [Kelvin]

        t0 (float or ndarray): Line temperature [Kelvin]

    Returns
        K2 (like input): How stimulated the emission is

    .. math::
        S(T) = S(T_0)K_1K_2 \\frac{Q(T_0)}{Q(T)}
    """
    return (1. - np.exp(- constants.planck * f0/(constants.boltzmann * t))) / \
        (1. - np.exp(- constants.planck * f0/(constants.boltzmann * t0)))
