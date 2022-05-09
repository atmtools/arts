import pyarts.arts as cxx
import test_functions as test

import numpy as np

x = cxx.ArrayOfIndex([1])
test.io(x, delete=True)
test.array(x)

x = cxx.ArrayOfIndex([1,2,3,4])
x = cxx.ArrayOfIndex(np.zeros(shape=(5), dtype=int))
assert np.all(np.array(x) == 0)

np.array(x, copy=False)[:] = 1
assert np.all(np.array(x) == 1)
