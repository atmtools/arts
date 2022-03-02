import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.ArrayOfPpath(1)
test.io(x, delete=True)
test.array(x)
