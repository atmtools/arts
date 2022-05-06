import pyarts.arts as cxx
import test_functions as test

x = cxx.ArrayOfSparse(1, cxx.Sparse())
test.io(x, delete=True)
test.array(x)
