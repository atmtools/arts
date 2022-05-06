import pyarts.arts as cxx
import test_functions as test

x = cxx.MapOfErrorCorrectedSuddenData()
test.io(x, delete=True)
