import pyarts.arts as cxx
import test_functions as test

import TestPropagationMatrix

x = cxx.ArrayOfPropagationMatrix([TestPropagationMatrix.x])
test.io(x, delete=True)
test.array(x)
