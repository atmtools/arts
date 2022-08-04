import numpy as np
import pyarts
import test_functions as test

ws = pyarts.workspace.Workspace()
ws.stokes_dim = 1
ws.f_grid = [1e9, 2e9, 3e9]

ws.starOff()
ws.starsAddSingleBlackbody(radius=20,
                       distance=2000,
                       temperature=5000,
                       latitude=10,
                       longitude=45)


star = ws.stars.value[0]

assert star.radius == 20
assert star.distance == 2000
assert star.latitude == 10
assert star.longitude == 45
assert np.isclose(star.spectrum[0, 0], 4.82602e-18, atol=1e-25)
assert np.isclose(star.spectrum[1, 0], 1.93040e-17, atol=1e-25)
assert np.isclose(star.spectrum[2, 0], 4.34338e-17, atol=1e-25)


x = pyarts.arts.ArrayOfStar(1, pyarts.arts.Star())
test.io(x, delete=True)
test.array(x)
