import pyarts.pyarts_cpp as cxx
import test_functions as test

import TestAgenda

x = cxx.ArrayOfAgenda([TestAgenda.x])

test.array(x)
