from pyarts.classes import (ArrayOfGriddedField2, ArrayOfString, GriddedField2, Matrix,
                            Vector, XsecRecord)

xr = XsecRecord()
xr2 = XsecRecord()

# VERSION 2
xr = XsecRecord()
xr.version = 2
xr.spec = "CFC11"
xr.fitminpressures = Vector([2, 3])
xr.fitmaxpressures = Vector([4, 5])
xr.fitmintemperatures = Vector([6, 7])
xr.fitmaxtemperatures = Vector([8, 9])
gf = GriddedField2()
gf.grids = [Vector([10, 11]), ArrayOfString(["p00", "p01", "p10", "p20"])]
gf.data = Matrix([[14, 15, 16, 17], [18, 19, 20, 21]])
xr.fitcoeffs = ArrayOfGriddedField2([gf, gf])

xr2.set(xr)
assert xr == xr2

xr3 = XsecRecord()
xr.savexml("tmp2.xr.xml", "ascii")
xr3.readxml("tmp2.xr.xml")
assert xr3 == xr
