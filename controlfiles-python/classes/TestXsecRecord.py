from pyarts.workspace import Workspace
from pyarts.classes.XsecRecord import XsecRecord
from pyarts.classes.Vector import Vector, ArrayOfVector
from pyarts.classes.Matrix import Matrix
from pyarts.classes.GriddedField2 import ArrayOfGriddedField2, GriddedField2
from pyarts.classes import from_workspace

xr = XsecRecord()
xr2 = XsecRecord()

# VERSION 2
xr = XsecRecord()
xr.version = 2
xr.spec = 1
xr.fitminpressures = Vector([2,3])
xr.fitmaxpressures = Vector([4,5])
xr.fitmintemperatures = Vector([6,7])
xr.fitmaxtemperatures = Vector([8,9])
gf = GriddedField2()
gf.grids=[Vector([10,11]),Vector([12,13])]
gf.data=Matrix([[14,15],[16,17]])
xr.fitcoeffs = ArrayOfGriddedField2([gf])

xr2.set(xr)
assert xr == xr2

xr3 = XsecRecord()
xr.savexml("tmp2.xr.xml", "ascii")
xr3.readxml("tmp2.xr.xml")
assert xr3 == xr
