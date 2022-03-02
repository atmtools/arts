import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.ArrayOfString(["OI"])
test.io(x, delete=False)
test.array(x)

x = cxx.ArrayOfString(["OI", cxx.String("AI")])
