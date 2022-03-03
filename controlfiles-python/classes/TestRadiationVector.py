import pyarts.pyarts_cpp as cxx
import test_functions as test

import numpy as np

x = cxx.RadiationVector(10, 4)
test.io(x, delete=True)

assert np.all(x[0] == 0)

x[0][:] = 2
assert np.all(x[0] == 2)
