import pyarts3 as pyarts
import numpy as np
import itertools

"""
Tests the internal wind jacobi with some limits considered

Note that the perturbations cannot be tested for all zeroes as
nor for both u and v zero. 
"""


ws = pyarts.Workspace()

dx = 1e-2

u = [10, 1, dx]
v = [10, 1, 0]
w = [10, 1, 0]
z = [30, 0, 90, 180]
a = [30, 0, 90, 170, 180]
f = [1e6, 1e9, 1e12]


COUNT = 0
for e in itertools.product(u, v, w, z, a, f):
    ws.atmospheric_point.wind = [e[0], e[1], e[2]]
    ws.ray_path_point.los = [e[3], e[4]]
    ws.frequency_grid = [e[5]]

    ws.frequency_gridWindShift()

    f0 = ws.frequency_grid[0] * 1.0
    df = 1.0 * ws.frequency_grid_wind_shift_jacobian * f0
    df_p = [0, 0, 0]

    for i in range(3):
        ws.ray_path_point.los = [e[3], e[4]]
        ws.atmospheric_point.wind = [e[0], e[1], e[2]]
        ws.atmospheric_point.wind[i] += dx
        ws.frequency_grid = [e[5]]
        ws.frequency_gridWindShift()
        f1 = ws.frequency_grid[0] * 1.0
        df_p[i] = (f1 - f0) / dx
    df_p = np.array(df_p)

    print(f"Wind: {[e[0], e[1], e[2]]}")
    print(f"LOS: {[e[3], e[4]]}")
    print(f"Frequency: {e[5]}")
    print(f"Jacobian: {df}")
    print(f"Perturbed: {df_p}")
    print(f"Test {COUNT}")
    assert np.allclose(df_p, df, rtol=1e-4) or np.allclose(
        df @ df, df_p @ df_p, rtol=1e-12
    )
    print()

    COUNT += 1
