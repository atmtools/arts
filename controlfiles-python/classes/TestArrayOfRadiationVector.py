import pyarts.arts as cxx
import test_functions as test

import TestRadiationVector

x = cxx.ArrayOfRadiationVector([TestRadiationVector.x])
test.io(x, delete=True)
test.array(x)
