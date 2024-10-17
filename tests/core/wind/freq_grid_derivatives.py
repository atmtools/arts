import pyarts
import numpy as np

ws = pyarts.Workspace()


winds = [
    [0, 0, 0],
    [10, 10, 10],
    # [10, 0, 0],
    # [0, 10, 0],
    [0, 0, 10],
    # [10, 10, 0],
    # [10, 0, 10],
    [0, 10, 10],
]

dx = 1e-3

zas = np.array([0, 45, 90, 180])
aas = np.array([0, 45, 90, 180, 189])
fs = np.array([1e6, 1e9, 1e12])

for wind in winds:
    print("wind:", wind)

    ws.atmospheric_point.temperature = 250
    ws.atmospheric_point.pressure = 1.0
    ws.atmospheric_point[pyarts.arts.SpeciesEnum.O2] = 0.2
    ws.atmospheric_point.wind = wind

    for za in zas:
        for aa in aas:
            for i in range(3):
                ws.frequency_grid = fs
                ws.ray_path_point.los = [za, aa]

                ws.atmospheric_point.wind = wind
                ws.frequency_gridWindShift()
                f0 = ws.frequency_grid * 1.0
                df0 = ws.frequency_grid_wind_shift_jacobian * 1.0

                ws.atmospheric_point.wind = wind
                ws.atmospheric_point.wind[i] += dx
                ws.frequency_grid = [1e6, 1e9, 1e12]
                ws.frequency_gridWindShift()
                f1 = ws.frequency_grid * 1.0
                df1 = ws.frequency_grid_wind_shift_jacobian * 1.0

                if i == 1:
                  print("uvw"[i], f"za: {za}, aa: {aa};", ((f1 - f0) / dx)[0] , ((f1 - f0) / dx)[0] / (fs * df0[i])[0], ((f1 - f0) / dx)[0] / (fs * df1[i])[0])
