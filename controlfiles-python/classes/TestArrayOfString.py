import pyarts.arts as cxx
import test_functions as test

x = cxx.ArrayOfString(["OI"])
test.io(x, delete=True)
test.array(x)

x = cxx.ArrayOfString(["OI", cxx.String("AI")])
