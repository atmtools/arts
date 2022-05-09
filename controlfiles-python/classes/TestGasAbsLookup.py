import pyarts.arts as cxx
import test_functions as test

x = cxx.GasAbsLookup()
test.io(x, delete=True)
