import pyarts.pyarts_cpp as cxx
import test_functions as test

import numpy as np

x = cxx.ArrayOfArrayOfVector([[[1,2,3]]])
test.io(x, delete=True)
test.array(x)
test.array_of_array(x)

x = cxx.ArrayOfArrayOfVector(np.zeros(shape=(3,3,3)))
assert np.all(np.array(x) == 0)

np.array(x[0][0])[:] = 1
assert np.all(np.array(x) == 0)

np.array(x[0][0], copy=False)[:] = 1
assert not np.all(np.array(x) == 0)
