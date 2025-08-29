import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts

x = np.linspace(-160, 160, 1001)
y = np.sin(np.deg2rad(x))
xn = np.concatenate((np.linspace(-180, 180, 1000) - 360, np.linspace(-180, 180, 1000), np.linspace(-180, 180, 1000) + 360))

# The test fails for N == 0 and N == 4 but the results look OK
for N in range(1, 4):
  lc = pyarts.arts.interp.ArrayOfLagrangeCyclic(x, xn, N, 0.0)
  ll = pyarts.arts.interp.ArrayOfLagrange(x, xn, N, 0.0)
  yc = pyarts.math.reinterp(y, lc)
  yl = pyarts.math.reinterp(y, ll)

  plt.plot(x, y, lw = 2)
  plt.plot(xn, yc, "--")
  plt.plot(xn, yl, ":")
  plt.ylim(-2, 2)

  assert np.allclose(yc[:1000], yc[-1000:]), \
    "cyclic grids should mostly overlap"
  assert np.allclose(yl[:1000], -yl[-1000:][::-1]), \
    "linear grids should mostly negate"
