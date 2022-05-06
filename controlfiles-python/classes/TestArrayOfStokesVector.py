import pyarts.arts as cxx
import test_functions as test

x = cxx.ArrayOfStokesVector(1, cxx.StokesVector())
test.io(x, delete=True)
test.array(x)
