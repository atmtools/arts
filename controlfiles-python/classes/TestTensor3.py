import pyarts.arts as cxx
import test_functions as test

import numpy as np

x = cxx.Tensor3([[[1, 2, 3]]])
test.io(x, delete=True)


x = cxx.Tensor3(np.zeros(shape=(3, 3, 3)))
assert np.all(np.array(x) == 0)

np.array(x)[:] = 1
assert np.all(np.array(x) == 0)

np.array(x, copy=False)[:] = 1
assert np.all(np.array(x) == 1)

x += 1
assert np.all(np.array(x) == 2)

x *= 2
assert np.all(np.array(x) == 4)
