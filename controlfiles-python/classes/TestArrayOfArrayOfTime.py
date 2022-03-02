import pyarts.pyarts_cpp as cxx
import test_functions as test


x = cxx.ArrayOfArrayOfTime(1, cxx.ArrayOfTime(1))
test.io(x, delete=True)
test.array(x)
test.array_of_array(x)
