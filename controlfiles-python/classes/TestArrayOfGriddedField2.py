import pyarts.pyarts_cpp as cxx
import test_functions as test

import TestGriddedField2

x = cxx.ArrayOfGriddedField2([TestGriddedField2.x])
test.io(x, delete=True)
test.array(x)