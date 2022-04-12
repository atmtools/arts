import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.SpeciesIsotopologueRatios()
test.io(x, delete=True)
