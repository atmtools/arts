"""GaussianRectangularAntenna samples a Gaussian beam on a rectangular grid of offsets

The samples are at all pairs of offsets towards increasing zenith and azimuth at
the boresight, weighted by the Gaussian gain of their angle from the boresight
times the solid angle of their cell.  So the weights sum to one, the samples
along the zenith axis stay in the vertical plane of the boresight, every sample
is at its angle from the boresight, and the beam has the second moment 2 std^2
of a Gaussian.  A vertical polarization response is not rotated along the
zenith axis, as that is the plane it is defined in.
"""

import numpy as np
import pyarts3 as pyarts

sensors = pyarts.arts.sensors

STD = 0.06  # [deg]
offsets = np.linspace(-5 * STD, 5 * STD, 201)
antenna = sensors.GaussianRectangularAntenna(offsets, offsets, STD)

pos = pyarts.arts.Vector3([800e3, 10.0, 20.0])
bore = np.array([130.0, 40.0])
ell = pyarts.arts.planets.Earth.ellipsoid

channel = sensors.DiracChannel(100e9)
channel.polarization.state = [1.0, 1.0, 0.0, 0.0]
sensor = antenna(channel, pos, bore, ell)

los = np.array(sensor.geometry)[:, 3:]
weights = np.array(sensor.elements[0].dense)[:, 0, :]
assert los.shape[0] == offsets.size**2

# The weights sum to one
assert np.isclose(weights[:, 0].sum(), 1.0, rtol=1e-14, atol=0)


def direction(l):
    za, aa = np.radians(l[..., 0]), np.radians(l[..., 1])
    return np.stack([np.sin(za) * np.sin(aa), np.sin(za) * np.cos(aa), np.cos(za)], axis=-1)


# The samples are by zenith offset and then azimuth offset; the middle azimuth offset is zero
z, a = np.meshgrid(offsets, offsets, indexing="ij")
z, a = z.ravel(), a.ravel()
axis = np.isclose(a, 0.0)
assert np.allclose(los[axis, 0], bore[0] + z[axis], rtol=0, atol=1e-9)
assert np.allclose(los[axis, 1], bore[1], rtol=0, atol=1e-9)

# Every sample is at its angle from the boresight
angle = np.degrees(np.arccos(np.clip(direction(los) @ direction(bore), -1, 1)))
expected = np.degrees(np.arctan(np.hypot(np.tan(np.radians(z)), np.tan(np.radians(a)))))
assert np.allclose(angle, expected, rtol=0, atol=1e-9)

# The second moment of the angle is that of a Gaussian of the angular offset, 2 std^2
second_moment = np.sum(weights[:, 0] * angle**2)
assert np.isclose(second_moment, 2 * STD**2, rtol=1e-4, atol=0), second_moment / (2 * STD**2) - 1

# A vertical response is not rotated in the vertical plane of the boresight
assert np.allclose(weights[axis, 1], weights[axis, 0], rtol=1e-12, atol=0)
