import pyarts.arts as cxx
import test_functions as test

x = cxx.ArrayOfScatteringMetaData(1, cxx.ScatteringMetaData())
test.io(x, delete=True)
test.array(x)
