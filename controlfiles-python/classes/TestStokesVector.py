import pyarts.pyarts_cpp as cxx
import test_functions as test

import numpy as np

x = cxx.StokesVector(50, 4)
test.io(x, delete=True)

assert np.all(np.array(x.data) == 0)

np.array(x.data)[:] = 1
assert np.all(np.array(x.data) == 0)

np.array(x.data, copy=False)[:] = 1
assert np.all(np.array(x.data) == 1)
