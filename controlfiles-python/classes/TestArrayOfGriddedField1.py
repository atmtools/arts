import pyarts.pyarts_cpp as cxx
import test_functions as test

import TestGriddedField1

x = cxx.ArrayOfGriddedField1([TestGriddedField1.x])
test.io(x, delete=True)
test.array(x)
