import pyarts.arts as cxx
import test_functions as test

import TestTransmissionMatrix

x = cxx.ArrayOfTransmissionMatrix([TestTransmissionMatrix.x])
test.io(x, delete=True)
test.array(x)
