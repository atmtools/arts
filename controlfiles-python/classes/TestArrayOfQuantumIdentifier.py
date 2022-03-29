import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.ArrayOfQuantumIdentifier(["H2O-161 J 1 1"])
test.io(x, delete=True)
test.array(x)
