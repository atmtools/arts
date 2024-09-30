# -*- coding: utf-8 -*-
"""
Implementation of functions related to sensor settings.

"""

import numpy as np


__all__ = ['get_f_backend_rel_width',
           'get_f_backend_const_width',
           ]


def get_f_backend_rel_width(f_start, f_end, bandwidth):
    """Compute backend frequencies with relative bandwidth.

    This function computes backend frequencies for a given frequency range and
    a relative bandwidth.

    Parameters:
        f_start (float):  beginning of frequency range [Hz]
        f_end (float): end of frequency range [Hz]
        bandwidth (float): relative bandwidth [dimensionless]

    Return:
        np.array, np.array: backend frequencies [Hz], channel widths [Hz]

    """
    if f_start <= 0:
        raise Exception('Start frequency must be > 0.')

    if f_start > f_end:
        raise Exception('End frequency has to be larger than start frequency.')

    f_backend = [f_start]
    while f_backend[-1] <= f_end:
        f_backend.append(f_backend[-1] * (bandwidth + 2) / (2 - bandwidth))

    # do not include last value in results as it exceeds f_end
    f_backend = np.array(f_backend[:-1])
    backend_bandwidth = f_backend * bandwidth

    return f_backend, backend_bandwidth


def get_f_backend_const_width(f_start, f_end, bandwidth):
    """Compute backend frequencies with constant bandwidth.

    This function computes backend frequencies for a given frequency range and
    a constant bandwidth.

    Parameters:
        f_start (float):  beginning of frequency range [Hz]
        f_end (float): end of frequency range [Hz]
        bandwidth (float): bandwidth [Hz]

    Return:
        np.array, np.array: backend frequencies [Hz], channel widths [Hz]

    """
    if f_start <= 0:
        raise Exception('Start frequency must be > 0.')

    if f_start > f_end:
        raise Exception('End frequency has to be larger than start frequency.')

    f_backend = [f_start]
    while f_backend[-1] <= f_end:
        f_backend.append(f_backend[-1] + bandwidth)

    # do not include last value in results as it exceeds f_end
    f_backend = np.array(f_backend[:-1])
    backend_bandwidth = np.array([bandwidth])

    return f_backend, backend_bandwidth
