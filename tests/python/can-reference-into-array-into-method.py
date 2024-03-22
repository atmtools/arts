"""
This test exsts because we are using shared pointers in ARTS to deal with the
workspace but pointers in the method interface to allow ws.method(a=b[0])-like
constructs.

It used to fail but should clearly work.
"""

import pyarts

ws = pyarts.Workspace()

ws.specs = pyarts.arts.ArrayOfArrayOfString([["O2"]])

ws.absorption_speciesSet(species=ws.specs[0])

x = pyarts.arts.ArrayOfArrayOfSpeciesTag()

ws.absorption_speciesSet(x, species=ws.specs[0])
