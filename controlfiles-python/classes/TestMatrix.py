import pyarts.pyarts_cpp as cxx
import test_functions as test

import numpy as np

x = cxx.Matrix([[1, 2, 3]])
test.io(x, delete=True)

y = cxx.Matrix([[1], [2], [3]])
assert (x @ y).shape == [1, 1]
assert (y @ x).shape == [3, 3]

z = cxx.Vector([1, 2, 3])
assert ((y @ x) @ z).shape == [3]

x = cxx.Matrix(np.zeros(shape=(3, 3)))
assert np.all(np.array(x) == 0)

np.array(x)[:] = 1
assert np.all(np.array(x) == 0)

np.array(x, copy=False)[:] = 1
assert np.all(np.array(x) == 1)

x += 1
assert np.all(np.array(x) == 2)

x *= 2
assert np.all(np.array(x) == 4)
