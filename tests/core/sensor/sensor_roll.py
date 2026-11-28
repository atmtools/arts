"""The roll of a SensorBuilder placement turns the whole instrument about its boresight

1. A roll turns the polarization response from v towards h, as
   Polarization.at turns it: a pencil beam with a V response rolled by chi
   receives [1, cos 2 chi, sin 2 chi, 0].
2. A roll of 0 is no roll, to the bit, for every kind of antenna pattern.
3. Looking at nadir, the frame of the instrument at roll 0 is that of the
   azimuth of the bore line of sight, so a placement at [180, aa] rolled by
   chi is the placement at [180, aa - chi] without a roll: the same lines of
   sight and the same weights, polarization included.  This holds for the
   pattern of rings, the rectangular pattern and the camera (whose pixels and
   beams turn), none with a sample exactly at nadir, where the azimuth of a
   line of sight is only a label.  A polarization that turned the other way
   from the pattern would fail it.
4. The number of rolls must be that of the placements.
"""

import numpy as np
import pyarts3 as pyarts

sensors = pyarts.arts.sensors
ell = pyarts.arts.Vector2([6378137.0, 6356752.314245])
pos = [pyarts.arts.Vector3([500e3, 10.0, 20.0])]

vertical = sensors.DiracChannel(100e9)
vertical.polarization.state = [1.0, 1.0, 0.0, 0.0]


def dense(sensor):
    return np.array([e.dense for e in sensor.elements])


# 1. The response turns from v towards h
chi = 30.0
pencil = pyarts.arts.SensorBuilder([vertical], sensors.PencilBeamAntenna())
rolled = pencil(pos, [pyarts.arts.Vector2([120.0, 40.0])], ell, roll=[chi])
weights = dense(rolled)[0, 0, 0]
expected = [1.0, np.cos(np.radians(2 * chi)), np.sin(np.radians(2 * chi)), 0.0]
assert np.allclose(weights, expected, rtol=0, atol=1e-14), (weights, expected)
turned = sensors.Polarization()
turned.state = [1.0, 1.0, 0.0, 0.0]
turned.scan_rotation = 1.0
assert np.allclose(weights, turned.at(chi), rtol=0, atol=1e-14), "the roll must turn as Polarization.at"

gain = np.array([[1.0, 0.7, 0.4, 0.2], [0.5, 0.3, 0.1, 0.05]])  # not symmetric about the boresight
antennas = {
    "pencil": sensors.PencilBeamAntenna(),
    "gaussian": sensors.GaussianAntenna([0.0, 0.5, 1.0], 0.4, azi_grid_size=6),
    "airy": sensors.GaussianAiryAntenna([0.0, 0.3, 0.6], 0.5, azi_grid_size=5),
    "gridded": sensors.GriddedAntennaPattern(
        sensors.AntennaPatternField(
            data=pyarts.arts.Matrix(gain),
            grids=(pyarts.arts.ZenGrid([0.5, 1.0]), pyarts.arts.AziGrid([0.0, 90.0, 180.0, 270.0])),
            grid_names=["zen", "azi"],
        )
    ),
    "rectangular": sensors.GaussianRectangularAntenna(
        zen_offsets=[-0.6, -0.2, 0.3, 0.7], azi_offsets=[-0.5, 0.1, 0.4], std=0.4
    ),
    "camera": sensors.CameraAntenna(
        sensors.PixelGrid(rows=[-0.8, 0.8], cols=[-0.5, 0.0, 0.5]),
        sensors.GaussianAntenna([0.0, 0.2], 0.1, azi_grid_size=3),
    ),
}

# 2. A roll of 0 is no roll
for kind, antenna in antennas.items():
    builder = pyarts.arts.SensorBuilder([vertical], antenna)
    los = [pyarts.arts.Vector2([130.0, 50.0])]
    plain, zero = builder(pos, los, ell), builder(pos, los, ell, roll=[0.0])
    assert np.array_equal(plain.geometry, zero.geometry), kind
    assert np.array_equal(dense(plain), dense(zero)), kind

# 3. At nadir a roll is a smaller azimuth of the bore line of sight
aa, chi = 40.0, 25.0
for kind in ["gridded", "rectangular", "camera"]:
    builder = pyarts.arts.SensorBuilder([vertical], antennas[kind])
    rolled = builder(pos, [pyarts.arts.Vector2([180.0, aa])], ell, roll=[chi])
    turned = builder(pos, [pyarts.arts.Vector2([180.0, aa - chi])], ell)
    los_r = np.asarray(rolled.geometry)[:, 3:]  # [alt, lat, lon, za, aa]
    los_t = np.asarray(turned.geometry)[:, 3:]
    assert np.all(los_r[:, 0] < 180.0) and np.all(los_t[:, 0] < 180.0), f"{kind}: a sample is at nadir"
    d_los = np.abs(np.remainder(los_r - los_t + 180.0, 360.0) - 180.0).max()
    d_weights = np.abs(dense(rolled) - dense(turned)).max() / np.abs(dense(turned)).max()
    print(f"{kind}: rolled by {chi} at [180, {aa}] vs [180, {aa - chi}]: max |dlos| {d_los:.1e} deg, "
          f"max |dweights| / max {d_weights:.1e}")
    assert d_los < 1e-9 and d_weights < 1e-9, kind

# 4. One roll per placement
try:
    pencil(pos * 2, [pyarts.arts.Vector2([120.0, 40.0])] * 2, ell, roll=[1.0])
    raise AssertionError("a roll per placement is required")
except RuntimeError:
    pass
