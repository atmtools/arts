import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.ArrayOfRetrievalQuantity(1, cxx.RetrievalQuantity())
# test.io(x, delete=True)
test.array(x)
