import pyarts.arts as cxx
import test_functions as test

import numpy as np

x = cxx.TransmissionMatrix(10, 4)
test.io(x, delete=True)

assert np.all(x[0] == np.diag([1,1,1,1]))

x[0][:] = 2
assert np.all(x[0] == 2)
