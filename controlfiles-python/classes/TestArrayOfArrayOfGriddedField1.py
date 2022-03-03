import pyarts.pyarts_cpp as cxx
import test_functions as test

import TestArrayOfGriddedField1

x = cxx.ArrayOfArrayOfGriddedField1([TestArrayOfGriddedField1.x])

test.io(x, delete=True)
test.array(x)
test.array_of_array(x)
