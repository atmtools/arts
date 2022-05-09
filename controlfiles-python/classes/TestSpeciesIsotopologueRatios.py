import pyarts.arts as cxx
import test_functions as test

x = cxx.SpeciesIsotopologueRatios()
test.io(x, delete=True)
