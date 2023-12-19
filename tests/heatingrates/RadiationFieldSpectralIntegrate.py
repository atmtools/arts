#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Manfred Brath

"""
import numpy as np

import pyarts

ws = pyarts.Workspace()

# Number of frequencies
N_f = 4

# Set test frequency grid
ws.f_grid = np.linspace(0, 2, N_f)  # Hz

# Set quadrature test weights
# Here we use linspace to make dummy quadrature weights
# The quadrature weights are normalized to the range of the frequency grid
quadrature_weights = np.linspace(1, np.pi, N_f)
quadrature_weights /= np.sum(quadrature_weights)
quadrature_weights *= ws.f_grid.value[-1] - ws.f_grid.value[0]

# Create  test spectral irradiance field
spectral_irradiance_field = np.zeros((N_f, 1, 1, 1, 2))
spectral_irradiance_field[:, 0, 0, 0, 0] = -np.linspace(0, 1, N_f)
spectral_irradiance_field[:, 0, 0, 0, 1] = np.linspace(0, 1, N_f)
ws.spectral_irradiance_field = spectral_irradiance_field

# Create test spectral radiance field
spectral_radiance_field = np.zeros((N_f, 1, 1, 1, 2, 1, 1))
spectral_radiance_field[:, 0, 0, 0, 0, 0, 0] = np.linspace(0, 1, N_f)
spectral_radiance_field[:, 0, 0, 0, 1, 0, 0] = np.linspace(0, 1, N_f)
ws.spectral_radiance_field = spectral_radiance_field

# Reference calculation
irradiance_field_ref = np.sum(
    spectral_irradiance_field
    * quadrature_weights[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis],
    0,
)

# For radiance field reference, we have to remove the last dimension of dummy_field_rad,
# because in the ARTS function the last dimension (stokes dimension) is also removed.
radiance_field_ref = np.sum(
    spectral_radiance_field[:, :, :, :, :, :, 0]
    * quadrature_weights[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis, np.newaxis],
    0,
)

# Do the integration for the spectral irradiance field
ws.RadiationFieldSpectralIntegrate(
    ws.irradiance_field, ws.f_grid, ws.spectral_irradiance_field, quadrature_weights
)

# Do the integration for the spectral radiance field
ws.RadiationFieldSpectralIntegrate(
    ws.radiance_field, ws.f_grid, ws.spectral_radiance_field, quadrature_weights
)


# Compare
assert np.allclose(
    ws.irradiance_field.value, irradiance_field_ref, 1e-6
), "spectral irradiance integration failed"
assert np.allclose(
    ws.radiance_field.value, radiance_field_ref, 1e-6
), "spectral radiance integration failed"
