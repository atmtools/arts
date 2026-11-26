import numpy as np
import pyarts3 as pyarts

# The subsurface is a function of the altitude relative to the local surface, so
# the same depth below an uneven surface gives the same subsurface values.

ws = pyarts.Workspace()
ws.surf_fieldEarth()
ws.surf_field["h"] = pyarts.arts.GeodeticField2(
    grids=[[-10.0, 10.0], [-20.0, 20.0]],
    grid_names=["latitude", "longitude"],
    data=[[100.0, 200.0], [300.0, 400.0]],
)

z = np.array([-3.0, -2.0, -1.0, 0.0])
ws.subsurf_field.bottom_depth = z[0]
ws.subsurf_field["t"] = pyarts.arts.GeodeticField3(
    grids=[z, [0.0], [0.0]],
    grid_names=["altitude", "latitude", "longitude"],
    data=(250.0 + 10.0 * z).reshape(len(z), 1, 1),
)

depths = [0.0, -0.5, -2.0]
for lat, lon in [(0.0, 0.0), (-10.0, -20.0), (5.0, 15.0)]:
    h = ws.surf_field["h"](lat, lon)
    ws.ray_point.pos = [h, lat, lon]
    ws.ray_point.los = [180.0, 0.0]
    ws.ray_pathFromPointAndDepth(depth_profile=depths)
    ws.subsurf_profileFromPath()
    t = [p.temperature for p in ws.subsurf_profile]
    assert np.allclose(t, 250.0 + 10.0 * np.array(depths), rtol=0, atol=1e-10), (h, t)
