"""measurement_sensorsAddSensorBuilder appends the sensor of a SensorBuilder

The sensor is the one the builder gives at the same position and bore line of
sight on the ellipsoid of surf_field.  It is appended to measurement_sensors and
named as the other measurement_sensorsAdd methods name theirs: by the given
name, or sensor-N for the N-th sensor, with its meta named the same.
"""

import numpy as np
import pyarts3 as pyarts

sensors = pyarts.arts.sensors

vertical = sensors.BoxChannel(pyarts.arts.AscendingGrid([89e9, 89.5e9, 90e9]))
vertical.polarization.state = [1.0, 1.0, 0.0, 0.0]
builder = pyarts.arts.SensorBuilder(
    [vertical, sensors.DiracChannel(183e9)], sensors.GaussianAntenna([0.0, 0.5, 1.0], 0.4, azi_grid_size=4)
)
pos = pyarts.arts.Vector3([800e3, 10.0, 20.0])
los = pyarts.arts.Vector2([150.0, 30.0])

ws = pyarts.Workspace()
ws.surf_fieldEarth()
ws.measurement_sensorsInit()
ws.measurement_sensorsAddSensorBuilder(pos=pos, los=los, sensor_builder=builder, name="radiometer")
ws.measurement_sensorsAddSensorBuilder(pos=pos, los=los, sensor_builder=builder)

expected = builder(pos, los, ws.surf_field.ellipsoid)
assert len(ws.measurement_sensors) == 2
for sensor, name in zip(ws.measurement_sensors, ["radiometer", "sensor-2"]):
    assert sensor.name == name, sensor.name
    assert sensor.meta.dataname == name, sensor.meta.dataname
    assert np.array_equal(sensor.freq_grid, expected.freq_grid)
    assert np.array_equal(sensor.geometry, expected.geometry)
    assert np.array_equal(sensor.meta.data, expected.meta.data)
    assert len(sensor.elements) == len(expected.elements)
    for got, want in zip(sensor.elements, expected.elements):
        assert np.array_equal(got.dense, want.dense)
