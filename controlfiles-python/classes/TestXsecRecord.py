from pyarts.workspace import Workspace
from pyarts.classes.XsecRecord import XsecRecord
from pyarts.classes.Vector import Vector, ArrayOfVector
from pyarts.classes import from_workspace

xr = XsecRecord()
xr2 = XsecRecord()

xr.spec = 1
xr.coeffs = 2
xr.ref_pressure = 3
xr.ref_temperature = 4
xr.fgrids = ArrayOfVector([Vector(5)])
xr.xsecs = ArrayOfVector([Vector(6)])
xr.temperature_slope = ArrayOfVector([Vector(7)])
xr.temperature_intersect = ArrayOfVector([Vector(8)])

xr2.set(xr)
assert xr == xr2
