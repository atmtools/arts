import pyarts.arts as cxx
import test_functions as test

x = cxx.ArrayOfArrayOfSpeciesTag(["H2O,H2O-MPM93"])

test.io(x, delete=True)
test.array(x)
test.array_of_array(x)

x = cxx.ArrayOfArrayOfSpeciesTag(["H2O", "H2O-PWR98"])
assert len(x) == 2
