import pyarts3 as pa
import numpy as np

NZA = 5
NAZ = 30
POS = [1, 0, 0]
LOS = [45, 30]
ELL = pa.arts.planets.Earth.ellipsoid

zen_grid = np.geomspace(1, 3, NZA) - 1
azi_grid = np.linspace(0, 360, NAZ + 1)[:-1]

ant1 = pa.arts.sensor.GaussianAiryAntenna(zen_grid, 30e-2, NAZ, "Iv")
ant2 = pa.arts.sensor.GaussianAiryAntenna(zen_grid + 1e-3, 30e-2, NAZ, "I")
ch = pa.arts.sensor.BoxChannel(1e9, 1e12, 1001)

obsel1 = ant1(ch, POS, LOS, ELL)
obsel2 = ant2(ch, POS, LOS, ELL)

poslos1 = np.asarray(obsel1.poslos)
poslos2 = np.asarray(obsel2.poslos)
weights1 = np.asarray(obsel1.weight_matrix.dense)
weights2 = np.asarray(obsel2.weight_matrix.dense)

pa.plots.plot(obsel1, point_spread=True)
pa.plots.plot(obsel2, point_spread=True, frame='local', type='contour', levels=10)

assert np.allclose(weights1.sum(), 2.0), \
    "Should be normalized to 2 (since polarized)"
assert np.allclose(weights2.sum(), 1.0), \
    "Should be normalized to 1 (since non-polarized)"
assert len(poslos1) == (NZA * NAZ - NAZ + 1), \
    "Should have one poslos per non-zero DZA, but only one for zero DZA - has zero DZA"
assert len(poslos2) == NZA * NAZ, \
    "Should have one poslos per non-zero DZA, but only one for zero DZA - has no zero DZA"

