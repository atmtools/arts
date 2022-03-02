import pyarts.pyarts_cpp as cxx
import test_functions as test

import TestArrayOfRadiationVector

x = cxx.ArrayOfArrayOfRadiationVector([TestArrayOfRadiationVector.x])

test.io(x, delete=True)
test.array(x)
test.array_of_array(x)
