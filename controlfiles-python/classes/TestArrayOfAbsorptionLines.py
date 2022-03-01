import pyarts.pyarts_cpp as cxx
import test_functions as test

import TestAbsorptionLines

x = cxx.ArrayOfAbsorptionLines([TestAbsorptionLines.x])

test.io(x, delete=True)
test.array(x)
