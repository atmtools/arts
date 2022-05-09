import pyarts.arts as cxx
import test_functions as test

x = cxx.ArrayOfXsecRecord(1, cxx.XsecRecord())
test.io(x, delete=True)
test.array(x)
