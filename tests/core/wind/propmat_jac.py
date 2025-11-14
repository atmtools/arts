#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 09:30:55 2024

@author: richard
"""

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

ws = pyarts.Workspace()

# Set up the frequency grid
line_f0 = 118750348044.712
f = np.linspace(-50e6, 50e6, 1000) + line_f0

# Remove the line frequencies += 1MHz because of numerical issues
f = f[np.where(abs(f - line_f0) > 1000e3)]
keys = "uvw"
wind = [10, 10, 10]
dx = 1e-1

# Load spectroscopic data
ws.abs_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.abs_bandsSelectFrequencyByLine(fmin=40e9, fmax=120e9)
ws.WignerInit()

# Standard surface and atmospheric setup
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atm_fieldIGRF(time="2000-03-11 14:39:37")
ws.atm_field["wind_u"] = wind[0]
ws.atm_field["wind_v"] = wind[1]
ws.atm_field["wind_w"] = wind[2]

# Set up a ray path point
ws.ray_point.los = [40, 20]
ws.ray_point.pos = [30e3, 0, 0]

for i in range(3):
    # Set up the wind field Jacobian
    ws.jac_targetsInit()
    ws.jac_targetsAddWindField(component=keys[i])
    ws.jac_targetsFinalize(measurement_sensor=[])

    # Reset
    ws.atm_point = ws.atm_field(*ws.ray_point.pos)

    # Original
    ws.freq_grid = f
    ws.freq_gridWindShift()
    ws.spectral_propmatInit()
    ws.spectral_propmatAddLines()
    ws.spectral_propmat_jacWindFix()

    # Keep analytical derivative and value
    x0 = ws.spectral_propmat * 1.0
    dd = ws.spectral_propmat_jac[0] * 1.0

    # Modify the wind field and calculate perturbed value
    ws.atm_point.wind[i] += dx
    ws.freq_grid = f
    ws.freq_gridWindShift()
    ws.spectral_propmatInit()
    ws.spectral_propmatAddLines()
    ws.spectral_propmat_jacWindFix()

    # Keep perturbed value and calculate perturbed derivative
    x1 = ws.spectral_propmat * 1.0
    d = (x1 - x0) / dx

    # Plot and compare the ratio between pertrubed and analytical derivatives
    plt.plot((f - line_f0) / 1e6, (d)[:, 0] / (dd)[:, 0], label=f"wind {keys[i]}-field")
    assert np.allclose(d[:, 0] / dd[:, 0], 1, rtol=1e-3)

plt.legend()
plt.xlabel("Frequency offset [MHz]")
plt.ylabel("Jacobian ratio perturbed/analytical [-]")
plt.title("Wind Jacobian ratio")
