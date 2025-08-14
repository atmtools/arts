import numpy as np
import pyarts

x = np.array([0.1, 0.4, 0.6, 0.9])
y = x

xn = 0.5
yn = pyarts.math.interp(y, pyarts.arts.interp.LagrangeCyclic(xn, x, 1))

assert np.isclose(0.5, yn)

x = np.array([-90, 90])
y = x

xn = 0.0
yn = pyarts.math.interp(y, pyarts.arts.interp.LagrangeCyclic(xn, x, 1))

assert np.isclose(0.0, yn)

xn = -180.0
yn = pyarts.math.interp(y, pyarts.arts.interp.LagrangeCyclic(xn, x, 1))

assert np.isclose(0.0, yn)
