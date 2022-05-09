import pyarts.arts as cxx
import test_functions as test

x = cxx.ArrayOfSingleScatteringData(1, cxx.SingleScatteringData())
# test.io(x, delete=True)
test.array(x)
