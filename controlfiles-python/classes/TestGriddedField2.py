import pyarts.arts as cxx
import test_functions as test

x = cxx.GriddedField2()
test.io(x, delete=True)
