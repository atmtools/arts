import pyarts.arts as cxx
import test_functions as test

import TestArrayOfStokesVector

x = cxx.ArrayOfArrayOfStokesVector([TestArrayOfStokesVector.x])

test.io(x, delete=True)
test.array(x)
test.array_of_array(x)
