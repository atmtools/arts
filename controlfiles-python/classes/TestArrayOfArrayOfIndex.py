import pyarts.pyarts_cpp as cxx
import test_functions as test

import numpy as np

x = cxx.ArrayOfArrayOfIndex([[1,2,3]])
test.io(x, delete=True)
test.array(x)
test.array_of_array(x)

x = cxx.ArrayOfArrayOfIndex([[1,2,3], [1,2]])
x = cxx.ArrayOfArrayOfIndex([np.array([1,2,3]), [1,2]])
x = cxx.ArrayOfArrayOfIndex(np.zeros(shape=(3, 6), dtype=int))
assert np.all(np.array(x) == 0)

np.array(x[0], copy=False)[:] = 1
assert not np.all(np.array(x) == 0)
