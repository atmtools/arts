import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

N = 100
f = np.logspace(-3, 20, N)

ws.measurement_sensorSimpleGaussian(
    freq_grid=f,
    std=1e3,
    pos=np.array([0, 1, 2]),
    los=np.array([3, 4]),
)

assert len(ws.measurement_sensor.unique_freq_grids()) == 1
assert len(ws.measurement_sensor.unique_poslos_grids()) == 1

ws.measurement_sensorAddSimpleGaussian(
    freq_grid=f,
    std=1e4,
    pos=np.array([0, 1, 2]),
    los=np.array([3, 4]),
)

assert len(ws.measurement_sensor.unique_freq_grids()) == 2
assert len(ws.measurement_sensor.unique_poslos_grids()) == 2

ws.measurement_sensor.collect_freq_grids()
assert len(ws.measurement_sensor.unique_freq_grids()) == 1
assert len(ws.measurement_sensor.unique_poslos_grids()) == 2

ws.measurement_sensor.collect_poslos_grids()
assert len(ws.measurement_sensor.unique_freq_grids()) == 1
assert len(ws.measurement_sensor.unique_poslos_grids()) == 1

ws.measurement_sensorAddSimpleGaussian(
    freq_grid=f,
    std=1e5,
    pos=np.array([0, 1, 2]),
    los=np.array([3, 4]),
)

ws.measurement_sensorAddSimpleGaussian(
    freq_grid=f + 1,
    std=1e5,
    pos=np.array([0, 1, 2]),
    los=np.array([3, 4]),
)

assert len(ws.measurement_sensor.unique_freq_grids()) == 3
assert len(ws.measurement_sensor.unique_poslos_grids()) == 3

ws.measurement_sensor.collect_freq_grids()
assert len(ws.measurement_sensor.unique_freq_grids()) == 1
assert len(ws.measurement_sensor.unique_poslos_grids()) == 3

fn = ws.measurement_sensor.savexml("sensor_count.xml")
ws.measurement_sensor.readxml("sensor_count.xml")

assert len(ws.measurement_sensor.unique_freq_grids()) == 1
assert len(ws.measurement_sensor.unique_poslos_grids()) == 3

ws.measurement_sensor.collect_poslos_grids()
assert len(ws.measurement_sensor.unique_poslos_grids()) == 1
assert len(ws.measurement_sensor.unique_freq_grids()) == 1

ws.measurement_sensor.collect_freq_grids()
assert len(ws.measurement_sensor.unique_poslos_grids()) == 1
assert len(ws.measurement_sensor.unique_freq_grids()) == 1
