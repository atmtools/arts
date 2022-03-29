import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.QuantumIdentifier("O2-66 v 0 0")
test.io(x, delete=True)
