import pyarts.arts as cxx
import test_functions as test

x = cxx.GriddedField6()
test.io(x, delete=True)
