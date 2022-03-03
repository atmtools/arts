import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.TelsemAtlas()
test.io(x, delete=True)
