import pyarts.pyarts_cpp as cxx
import test_functions as test

import TestArrayOfPropagationMatrix

x = cxx.ArrayOfArrayOfPropagationMatrix([TestArrayOfPropagationMatrix.x])

test.io(x, delete=True)
test.array(x)
test.array_of_array(x)
