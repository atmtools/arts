import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.ArrayOfStokesVector(1)
test.io(x, delete=False)
test.array(x)
