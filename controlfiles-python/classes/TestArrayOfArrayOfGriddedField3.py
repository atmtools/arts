import pyarts.pyarts_cpp as cxx
import test_functions as test

import TestGriddedField3

x = cxx.ArrayOfArrayOfGriddedField3([TestGriddedField3.x])

test.io(x, delete=True)
test.array(x)
test.array_of_array(x)
