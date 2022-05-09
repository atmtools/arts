import pyarts.arts as cxx
import test_functions as test

import TestGriddedField3

x = cxx.ArrayOfGriddedField3([TestGriddedField3.x])
test.io(x, delete=True)
test.array(x)
