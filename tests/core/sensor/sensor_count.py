import pyarts3 as pyarts
import numpy as np

ws = pyarts.Workspace()

N = 100
f = np.logspace(-3, 20, N)
pos = np.array([0, 1, 2])
los = np.array([3, 4])

ws.measurement_sensorsSimpleGaussian(freq_grid=f, std=1e3, pos=pos, los=los)

# Each Add appends one sensor, with one element per frequency, named after its position unless a name is given
assert len(ws.measurement_sensors) == 1
assert ws.measurement_sensors[0].name == "sensor-1"
assert len(ws.measurement_sensors[0].elements) == N

ws.measurement_sensorsAddSimpleGaussian(freq_grid=f, std=1e4, pos=pos, los=los, name="wide")
ws.measurement_sensorsAddSimpleGaussian(freq_grid=f + 1, std=1e5, pos=pos, los=los)

assert [s.name for s in ws.measurement_sensors] == ["sensor-1", "wide", "sensor-3"]
assert [s.meta.dataname for s in ws.measurement_sensors] == ["sensor-1", "wide", "sensor-3"]
assert [len(s.geometry) for s in ws.measurement_sensors] == [1, 1, 1]

expected = pyarts.arts.ArrayOfSensor(ws.measurement_sensors)
ws.measurement_sensors.savexml("sensor_count.xml")
ws.measurement_sensors.readxml("sensor_count.xml")

assert [s.name for s in ws.measurement_sensors] == [s.name for s in expected]
for got, exp in zip(ws.measurement_sensors, expected):
    assert np.array_equal(got.freq_grid, exp.freq_grid)
    assert np.array_equal(got.geometry, exp.geometry)
    assert len(got.elements) == len(exp.elements)
    assert np.array_equal(got.elements[N // 2].dense, exp.elements[N // 2].dense)

# All three sensors look along the same line of sight, so each makes its own agenda call for it
report = ws.measurement_sensors.report()
print(report)
assert "wide" in report
assert "Geometries simulated by more than one sensor: 1, costing 2 agenda calls" in report

# The part of a measurement vector of each sensor, with the structure of its meta
y = np.arange(3.0 * N)
parts = ws.measurement_sensors.split(y)
assert [len(part) for part in parts] == [N, N, N]
assert np.array_equal(parts[1], y[N:2 * N])
assert np.array_equal(ws.measurement_sensors[1].reshape(parts[1]).data, y[N:2 * N])
assert np.array_equal(ws.measurement_sensors[1].reshape(parts[1]).grids[0], f)
