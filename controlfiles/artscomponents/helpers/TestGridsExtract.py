#
# Tests extraction of (computational) atmospheric grids from data grids. Also
# crops this grids to a user-defined region and tests functionality by
# interpolating the considered atmospheric fields on the extracted grids.
#
# Note: Cropping does not (yet?) take into account cyclicity of longitudinal
# grid (i.e. we can only crop within the data grid, which here is 0-360deg. That
# is, a cropped region, e.g., of -30..30 actually covering 330..360=0..30 is not
# possible (yet?), but -30..30 setting will result in 0..30deg grid.
#
# Jana Mendrok 2013-10-04

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
# create the bunch of local variables we are going to use later on
ws.NumericCreate("pmin")
ws.NumericCreate("pmax")
ws.NumericCreate("latmin")
ws.NumericCreate("latmax")
ws.NumericCreate("lonmin")
ws.NumericCreate("lonmax")
ws.StringCreate("caseext")
ws.StringCreate("casename")
ws.StringCreate("atmcase")
ws.StringCreate("Bname")
ws.GriddedField3Create("grid_field")
ws.NumericSet(ws.pmin, 10.0)
ws.NumericSet(ws.pmax, 200000.0)
ws.NumericSet(ws.latmin, -90.0)
ws.NumericSet(ws.latmax, -63.0)
ws.NumericSet(ws.lonmin, -17.0)
ws.NumericSet(ws.lonmax, 177.0)
ws.AtmosphereSet3D()
# set basic case folder
ws.StringCreate("basename")
ws.StringSet(ws.basename, "planets/Jupiter/MPS/")
# atmcase name
ws.StringSet(ws.atmcase, "Jupiter.mean")
ws.abs_speciesSet(species=["CH4", "H2", "H2-CIA-H2-0", "CH4-CIA-CH4-0"])
# construct atmcase name
ws.Append(ws.basename, ws.atmcase)
ws.StringSet(ws.caseext, "/")
ws.Append(ws.basename, ws.caseext)
ws.Append(ws.basename, ws.atmcase)
ws.AtmRawRead(basename=ws.basename)
ws.StringSet(ws.Bname, "planets/Jupiter/Khurana/Khurana.B_w.xml.gz")
ws.ReadXML(ws.grid_field, ws.Bname)
# now derive p_grid from given raw z_field and regrid atm fields to this
ws.p_gridFromZRaw(no_negZ=0)
ws.lat_gridFromRawField(field_raw=ws.grid_field)
ws.lon_gridFromRawField(field_raw=ws.grid_field)
ws.VectorCrop(ws.p_grid, ws.p_grid, ws.pmin, ws.pmax)
ws.VectorCrop(ws.lat_grid, ws.lat_grid, ws.latmin, ws.latmax)
ws.VectorCrop(ws.lon_grid, ws.lon_grid, ws.lonmin, ws.lonmax)
# Print( p_grid )
# Print( lat_grid )
# Print( lon_grid )
ws.AtmFieldsCalcExpand1D(vmr_zeropadding=1)
