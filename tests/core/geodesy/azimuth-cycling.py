"""East/west and north/south decisions do not depend on the azimuth range.

A line-of-sight azimuth may be given in [-180, 180] or in AziGrid's
[0, 360).  The geodesy routines cycle it to [0, 360) before deciding whether
a ray points east (0, 180) or west (180, 360), so the same direction in either
range must give the same answer.

The reference for the longitude crossing is a closed form: a horizontal ray
at the equator heading east from longitude 0 is the line x = a + h,
y = l, so it crosses the meridian of longitude L at l = (a + h) tan(L).
"""

import numpy as np
import pyarts3 as pyarts

geo = pyarts.arts.geodetic
ell = [6378137.0, 6356752.314245]
h = 1000.0
pos = [h, 0.0, 0.0]


def crossings(los):
    ecef, decef = geo.geodetic_los2ecef(pos, los, ell)
    return (
        geo.intersection_longitude(ecef, decef, pos, los, 1.0),
        geo.intersection_longitude(ecef, decef, pos, los, -1.0),
        geo.intersection_latitude(ecef, decef, pos, los, ell, 1.0),
        geo.intersection_latitude(ecef, decef, pos, los, ell, -1.0),
    )


east = (ell[0] + h) * np.tan(np.radians(1.0))
for aa in (90.0, 450.0, -270.0):
    lon_east, lon_west, *_ = crossings([90.0, aa])
    assert abs(lon_east - east) < 1e-12 * east, (aa, lon_east, east)
    assert lon_west == -1, (aa, lon_west)

for aa in (270.0, -90.0):
    lon_east, lon_west, *_ = crossings([90.0, aa])
    assert lon_east == -1, (aa, lon_east)
    assert abs(lon_west - east) < 1e-12 * east, (aa, lon_west, east)

# The same direction in both ranges, in every quadrant and on the axes: the
# same crossings (-1 for none), to the round-off of sin and cos of the azimuth
for aa in (0.0, 30.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0):
    for za in (60.0, 90.0, 120.0):
        reference = np.array(crossings([za, aa]))
        for other in (aa - 360.0, aa + 360.0):
            got = np.array(crossings([za, other]))
            assert np.array_equal(got == -1, reference == -1), (za, aa, other, got, reference)
            assert np.allclose(got, reference, rtol=1e-12, atol=0), (za, aa, other, got, reference)
