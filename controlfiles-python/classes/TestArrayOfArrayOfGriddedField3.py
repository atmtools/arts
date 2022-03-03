import pyarts.pyarts_cpp as cxx
import test_functions as test

import TestArrayOfGriddedField3

x = cxx.ArrayOfArrayOfGriddedField3([TestArrayOfGriddedField3.x])

test.io(x, delete=True)
test.array(x)
test.array_of_array(x)
