import pyarts.arts as cxx
import test_functions as test

import numpy as np

x = cxx.ArrayOfTensor4([[[[[1,2,3]]]]])
test.io(x, delete=True)
test.array(x)

x = cxx.ArrayOfTensor4(np.zeros(shape=(3, 3, 3, 3, 3)))
assert np.all(np.array(x) == 0)

np.array(x[0])[:] = 1
assert np.all(np.array(x) == 0)

np.array(x[0], copy=False)[:] = 1
assert not np.all(np.array(x) == 0)
