import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.ArrayOfSpeciesTag("H2O")
x = cxx.ArrayOfSpeciesTag(["H2O"])
test.io(x, delete=True)
test.array(x)
