import numpy as np
import pyarts3 as pyarts
import time

start_time = time.time()

NCHANNELS = 10_000
IF_LOW = 5.5e9
IF_HIGH = 6.5e9
LO = 600e9
APERTURE_DIAMETER = 30e-2
POS = [600e3, 0.0, 0.0]
LOS = [45.0, 0.0]
ZEN = np.linspace(0.0, 0.3, 7)
AZI = np.linspace(0.0, 360.0, 9)[:-1]


if_offsets = pyarts.arts.AscendingGrid(np.linspace(IF_LOW, IF_HIGH, NCHANNELS))
print(f"if_offsets:   {time.time() - start_time:.2f} seconds")

spectrometer = pyarts.arts.sensor.Spectrometer(
    pyarts.arts.sensor.DiracChannel(),
    if_offsets,
)
print(f"spectrometer: {time.time() - start_time:.2f} seconds")

backend = pyarts.arts.sensor.HeterodyneFrequencyRange(
    LO,
    np.array([LO + IF_LOW, LO + IF_HIGH]),
)
print(f"backend:      {time.time() - start_time:.2f} seconds")

antenna = pyarts.arts.sensor.GaussianAiryAntenna(
    ZEN,
    APERTURE_DIAMETER,
    8,
    "I",
)
print(f"antenna:      {time.time() - start_time:.2f} seconds")

sensor = pyarts.arts.sensor.Builder(antenna, spectrometer, backend)
print(f"sensor:       {time.time() - start_time:.2f} seconds")

ws = pyarts.Workspace()

ws.measurement_sensor, ws.measurement_sensor_meta = sensor(POS, LOS, pyarts.arts.planets.Earth.ellipsoid)

assert len(ws.measurement_sensor.unique_freq_grids()) == 1
assert len(ws.measurement_sensor[0].f_grid) == NCHANNELS

np.testing.assert_allclose(
    np.asarray(ws.measurement_sensor[0].f_grid)[[0, -1]],
    np.array([LO + IF_LOW, LO + IF_HIGH]),
    atol=0.0,
    rtol=0.0,
)
