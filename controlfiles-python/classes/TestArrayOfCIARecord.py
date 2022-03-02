import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.ArrayOfCIARecord()
test.io(x, delete=True)
