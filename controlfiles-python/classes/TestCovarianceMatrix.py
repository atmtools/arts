import pyarts.arts as cxx
import test_functions as test

x = cxx.CovarianceMatrix()
test.io(x, delete=True)
