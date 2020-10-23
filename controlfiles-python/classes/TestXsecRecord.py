from pyarts.workspace import Workspace
from pyarts.classes.XsecRecord import XsecRecord
from pyarts.classes.Vector import Vector, ArrayOfVector
from pyarts.classes.Matrix import Matrix
from pyarts.classes.GriddedField2 import ArrayOfGriddedField2, GriddedField2
from pyarts.classes import from_workspace

xr = XsecRecord()
xr2 = XsecRecord()

# VERSION 1
xr.version = 1
xr.spec = 1
xr.coeffs = 2
xr.ref_pressure = 3
xr.ref_temperature = 4
xr.fgrids = ArrayOfVector([Vector([5,6])])
xr.xsecs = ArrayOfVector([Vector([7,8])])
xr.temperature_slope = ArrayOfVector([Vector([9,10])])
xr.temperature_intersect = ArrayOfVector([Vector([11,12])])
xr.fitminpressures = Vector([2,3])

xr2.set(xr)
assert xr == xr2

xr3 = XsecRecord()
xr.savexml("tmp1.xr.xml", "ascii")
xr3.readxml("tmp1.xr.xml")


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
