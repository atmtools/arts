import pyarts.pyarts_cpp as cxx
import test_functions as test

import TestArrayOfTransmissionMatrix

x = cxx.ArrayOfArrayOfTransmissionMatrix([TestArrayOfTransmissionMatrix.x])

test.io(x, delete=True)
test.array(x)
test.array_of_array(x)
