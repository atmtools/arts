#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 09:30:55 2024

@author: richard
"""

import pyarts
import numpy as np

import matplotlib.pyplot as plt

ws = pyarts.Workspace()

line_f0 = 118750348044.712
f = np.linspace(-25e6, 25e6, 1001) + line_f0

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequency(fmin=40e9, fmax=120e9, by_line=1)
ws.WignerInit()

ws.surface_field
ws.atmospheric_field["wind_w"] = 0
ws.jacobian_targetsInit()
ws.jacobian_targetsAddWindField(component="w")
ws.jacobian_targetsFinalize()

wind = [10, 10, 10]
dx = 0.1

ws.atmospheric_point.temperature = 250
ws.atmospheric_point.pressure = 1.0
ws.atmospheric_point[pyarts.arts.SpeciesEnum.O2] = 0.2
ws.atmospheric_point.wind = wind

ws.ray_path_point.los = [40, 0]
ws.ray_path_point.pos = [90e3, 0, 0]

ws.frequency_grid = pyarts.arts.frequency_shift(f, ws.ray_path_point, ws.atmospheric_point)
print(ws.frequency_grid)
ws.propagation_matrixInit()
ws.frequency_gridWindShift()
ws.propagation_matrixAddLines()
print(ws.propagation_matrix_jacobian)
ws.propagation_matrix_jacobianWindFix()
print(ws.propagation_matrix_jacobian)

x0 = ws.propagation_matrix * 1.0
dd = ws.propagation_matrix_jacobian[0] * 1.0

ws.atmospheric_point.wind[2] += dx

ws.frequency_grid = pyarts.arts.frequency_shift(f, ws.ray_path_point, ws.atmospheric_point)
print(ws.frequency_grid)
ws.propagation_matrixInit()
ws.propagation_matrixAddLines()
ws.propagation_matrix_jacobianWindFix()

x1 = ws.propagation_matrix * 1.0
d = (x1 - x0) / dx


plt.semilogy((f - line_f0)/1e6, abs(d)[:, 0])
plt.semilogy((f - line_f0)/1e6, abs(dd)[:, 0])
plt.show()
plt.plot((f - line_f0)/1e6, (d)[:, 0] / (dd)[:, 0])
plt.show()
assert np.allclose(d, dd, rtol = 1e-3)
