import pyarts3 as pa
import numpy as np

NZA = 5
NAZ = 30
POS = [1, 0, 0]
LOS = [45, 30]
ELL = pa.arts.planets.Earth.ellipsoid

zen_grid = np.geomspace(1, 3, NZA) - 1
azi_grid = np.linspace(0, 360, NAZ + 1)[:-1]

ant1 = pa.arts.sensors.GaussianAiryAntenna(zen_grid, 30e-2, NAZ)
ant2 = pa.arts.sensors.GaussianAiryAntenna(zen_grid + 1e-3, 30e-2, NAZ)
ch1 = pa.arts.sensors.BoxChannel(1e9, 1e12, 1001)
ch1.polarization = "Iv"
ch2 = pa.arts.sensors.BoxChannel(1e9, 1e12, 1001)

sensor1 = ant1(ch1, POS, LOS, ELL)
sensor2 = ant2(ch2, POS, LOS, ELL)

poslos1 = np.asarray(sensor1.geometry)
poslos2 = np.asarray(sensor2.geometry)
weights1 = np.asarray(sensor1.elements[0].dense)
weights2 = np.asarray(sensor2.elements[0].dense)

pa.plots.plot(sensor1, point_spread=True)
pa.plots.plot(sensor2, point_spread=True, frame='local', type='contour', levels=10)

assert np.allclose(weights1[..., 0].sum(), 1.0), \
    "Should receive unpolarized radiation with weight 1"
assert np.allclose(weights2[..., 0].sum(), 1.0), \
    "Should receive unpolarized radiation with weight 1"
assert np.allclose(weights2[..., 1:], 0.0), \
    "An unpolarized channel has no polarized weights"
assert 0.99 < weights1[..., 1].sum() <= 1.0, \
    "A vertical channel off nadir receives almost only vertical polarization in its lines of sight"
assert len(poslos1) == (NZA * NAZ - NAZ + 1), \
    "Should have one poslos per non-zero DZA, but only one for zero DZA - has zero DZA"
assert len(poslos2) == NZA * NAZ, \
    "Should have one poslos per non-zero DZA, but only one for zero DZA - has no zero DZA"

