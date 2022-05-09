import pyarts.arts as cxx
import test_functions as test

import TestGriddedField4

x = cxx.ArrayOfGriddedField4([TestGriddedField4.x])
test.io(x, delete=True)
test.array(x)
