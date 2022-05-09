import pyarts.arts as cxx
import test_functions as test

import numpy as np

x = cxx.Vector([1, 2, 3])
test.io(x, delete=True)


x = cxx.Vector(np.zeros(shape=(3)))
assert np.all(np.array(x) == 0)

np.array(x)[:] = 1
assert np.all(np.array(x) == 0)

np.array(x, copy=False)[:] = 1
assert np.all(np.array(x) == 1)

x += 1
assert np.all(np.array(x) == 2)

x *= 2
assert np.all(np.array(x) == 4)

assert x @ x == 48
