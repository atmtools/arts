import numpy as np
from pyarts.workspace import Workspace
from pyarts.classes.Star import Star
from pyarts.classes import from_workspace

ws = Workspace()
ws.stokes_dim = 1
f_grid = from_workspace(ws.f_grid)
f_grid.data = [1e9, 2e9, 3e9]
stars = from_workspace(ws.stars)
stars.data = []

ws.starBlackbodySimple(radius=20,
                       distance=2000,
                       temperature=5000,
                       latitude=10,
                       longitude=45)

# Reach for the stars
stars = from_workspace(ws.stars)
star = stars[0]
star2 = Star()
star2.set(star)

assert star == star2

assert star.radius == 20
assert star.distance == 2000
assert star.latitude == 10
assert star.longitude == 45
assert np.isclose(star.spectrum.data[0, 0], 4.82602e-18, atol=1e-25)
assert np.isclose(star.spectrum.data[1, 0], 1.93040e-17, atol=1e-25)
assert np.isclose(star.spectrum.data[2, 0], 4.34338e-17, atol=1e-25)

star3 = Star()
star.savexml("tmp.star.xml", "binary")
star3.readxml("tmp.star.xml")
assert star == star3
