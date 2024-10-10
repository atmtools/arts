#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 09:30:55 2024

@author: richard
"""

import pyarts
import numpy as np

ws = pyarts.Workspace()

line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(5e6, 6e6, 3) + line_f0

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequency(fmin=40e9, fmax=120e9, by_line=1)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

B = [-3.132846e-06, 2.62680294e-05, 1.39844339e-05]
dx = 1e-11

B1 = [b for b in B]
B1[2] += dx

ws.surface_field
ws.atmospheric_field["mag_w"] = B[2]
ws.jacobian_targetsInit()
ws.jacobian_targetsAddMagneticField(component="w")
ws.jacobian_targetsFinalize()

ws.atmospheric_point.temperature = 250
ws.atmospheric_point.pressure = 1.0
ws.atmospheric_point[pyarts.arts.SpeciesEnum.O2] = 0.2
ws.atmospheric_point.mag = B

ws.ray_path_point.los = [40, 0]
ws.ray_path_point.pos = [90e3, 0, 0]

ws.propagation_matrixInit()
ws.propagation_matrixAddLines()

x0 = ws.propagation_matrix * 1.0
dd = ws.propagation_matrix_jacobian[0] * 1.0

ws.propagation_matrixInit()
ws.propagation_matrixAddLines()

ws.atmospheric_point.mag = B1

ws.propagation_matrixInit()
ws.propagation_matrixAddLines()

x1 = ws.propagation_matrix * 1.0
d = (x1 - x0) / dx

assert np.allclose(d, dd, rtol = 1e-3)
