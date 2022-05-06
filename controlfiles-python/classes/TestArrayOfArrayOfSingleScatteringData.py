import pyarts.arts as cxx
import test_functions as test

x = cxx.ArrayOfArrayOfSingleScatteringData()
test.io(x, delete=True)
