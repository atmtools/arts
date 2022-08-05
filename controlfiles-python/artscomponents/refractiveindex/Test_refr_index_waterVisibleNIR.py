#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Manfred Brath
"""

import numpy as np
import pyarts as pa


# =============================================================================
# %% define some useful conversion functions
# =============================================================================

def wavelength_in_microns2frequency(wavelength):
    '''
    Converts wavelength in µm to frequency in Hz.

    Args:
        wavelength (ndarray): Wavelength in µm.

    Returns:
        ndarry: Frequency in Hz.

    '''

    # speed of light
    c = pa.arts.constant.c  # [m/s]

    return c / (wavelength * 1e-6)  # [Hz]


# =============================================================================
# %% open workspace
# =============================================================================


ws = pa.workspace.Workspace()
ws.verbositySetScreen(level=2)

# Define wavelength in µm
wavelength = np.array([0.488])

data_f_grid = wavelength_in_microns2frequency(wavelength)

# Define temperature
data_t_grid = np.array([283.15])

# Set density
density_water = np.array([1000.])

ws.Touch(ws.complex_refr_index)
ws.refr_index_waterVisibleNIR(data_f_grid=data_f_grid,
                              data_t_grid=data_t_grid,
                              density_water=density_water)

complex_refr_index = ws.complex_refr_index.value

# Reference data
complex_refr_index_REFERENCE = pa.arts.GriddedField3([[614328807377049], [283.15], ["real", "imaginary"]],
                                                     [[[1.33821937010893, 0]]], ["Frequency", "Temperature", "Complex"])
# Compare with reference
ws.Compare(complex_refr_index, complex_refr_index_REFERENCE, 1e-6)
