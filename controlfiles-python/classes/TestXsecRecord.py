from pyarts.workspace import Workspace
from pyarts.classes.XsecRecord import XsecRecord
from pyarts.classes.Vector import Vector, ArrayOfVector
from pyarts.classes import from_workspace

xr = XsecRecord()
xr2 = XsecRecord()

xr.spec = "H2O"
xr.coeffs = 2
xr.ref_pressure = 3
xr.ref_temperature = 4
xr.fgrids = ArrayOfVector([Vector([5,6])])
xr.xsecs = ArrayOfVector([Vector([7,8])])
xr.temperature_slope = ArrayOfVector([Vector([9,10])])
xr.temperature_intersect = ArrayOfVector([Vector([11,12])])

xr2.set(xr)
assert xr == xr2

xr3 = XsecRecord()
xr.savexml("tmp.xr.xml", "binary")
xr3.readxml("tmp.xr.xml")
assert xr3 == xr
