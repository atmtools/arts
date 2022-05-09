import pyarts.arts as cxx
import test_functions as test

x = cxx.ArrayOfArrayOfString([["OI"]])
test.io(x, delete=True)
test.array(x)
test.array_of_array(x)
