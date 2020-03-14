import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.GriddedField3Create("gf_ref")
ws.GriddedField3Create("gf_regridded")
ws.GriddedField3Create("gf")
ws.StringCreate("fname")
ws.NumericCreate("maxabsdiff")
ws.NumericSet(ws.maxabsdiff, 1e-08)
ws.ReadXML(out=ws.gf, filename="gf_data.xml")
########## GriddedFieldPRegrid ##########
#
# New grid is inside the old grid
#
ws.StringSet(ws.fname, "gfREFERENCE_0p_none.xml")
ws.VectorNLogSpace(ws.p_grid, 20, 6000.0, 2000.0)
ws.GriddedFieldPRegrid(ws.gf_regridded, ws.p_grid, ws.gf, 1, 0)
# WriteXML("ascii", gf_regridded, fname)
ws.ReadXML(ws.gf_ref, ws.fname)
ws.Compare(ws.gf_regridded, ws.gf_ref, ws.maxabsdiff)
#
# New grid is larger than the old grid
#
ws.StringSet(ws.fname, "gfREFERENCE_0p_both_sides.xml")
ws.VectorNLogSpace(ws.p_grid, 20, 15000.0, 100.0)
ws.GriddedFieldPRegrid(ws.gf_regridded, ws.p_grid, ws.gf, 1, 1)
# WriteXML("ascii", gf_regridded, fname)
ws.ReadXML(ws.gf_ref, ws.fname)
ws.Compare(ws.gf_regridded, ws.gf_ref, ws.maxabsdiff)
#
# Minimum pressure in new grid is lower than in the old grid
#
ws.StringSet(ws.fname, "gfREFERENCE_0p_bottom.xml")
ws.VectorNLogSpace(ws.p_grid, 20, 5000.0, 100.0)
ws.GriddedFieldPRegrid(ws.gf_regridded, ws.p_grid, ws.gf, 1, 1)
# WriteXML("ascii", gf_regridded, fname)
ws.ReadXML(ws.gf_ref, ws.fname)
ws.Compare(ws.gf_regridded, ws.gf_ref, ws.maxabsdiff)
#
# Maximum pressure in the new grid is higher than in the old grid
#
ws.StringSet(ws.fname, "gfREFERENCE_0p_top.xml")
ws.VectorNLogSpace(ws.p_grid, 20, 15000.0, 5000.0)
ws.GriddedFieldPRegrid(ws.gf_regridded, ws.p_grid, ws.gf, 1, 1)
# WriteXML("ascii", gf_regridded, fname)
ws.ReadXML(ws.gf_ref, ws.fname)
ws.Compare(ws.gf_regridded, ws.gf_ref, ws.maxabsdiff)
########## GriddedFieldZToPRegrid ##########
ws.ReadXML(ws.gf, "gf.xml")
ws.GriddedField3Create("gf_z")
ws.ReadXML(out=ws.gf_z, filename="gf_data_z.xml")
ws.VectorSet(ws.lat_grid, np.array([0.0]))
ws.VectorSet(ws.lon_grid, np.array([0.0]))
ws.GriddedField3Create("zraw_regridded")
#
# New grid is inside the old grid
#
ws.StringSet(ws.fname, "gfREFERENCE_0p_none.xml")
ws.VectorNLogSpace(ws.p_grid, 20, 6000.0, 2000.0)
ws.GriddedFieldPRegrid(ws.zraw_regridded, ws.p_grid, ws.gf, 1, 0)
ws.FieldFromGriddedField(
    ws.z_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.zraw_regridded
)
ws.GriddedFieldZToPRegrid(
    ws.gf_regridded, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.z_field, ws.gf_z, 1, 1
)
ws.ReadXML(ws.gf_ref, ws.fname)
ws.Compare(ws.gf_regridded, ws.gf_ref, ws.maxabsdiff)
#
# New grid is larger than the old grid
#
ws.StringSet(ws.fname, "gfREFERENCE_z_0p_both_sides.xml")
ws.ReadXML(ws.p_grid, "p_grid_both.xml")
ws.ReadXML(ws.z_field_raw, "z_both.xml")
ws.FieldFromGriddedField(
    ws.z_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.z_field_raw
)
ws.GriddedFieldZToPRegrid(
    ws.gf_regridded, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.z_field, ws.gf_z, 1, 1
)
# WriteXML("ascii", gf_regridded, fname)
ws.ReadXML(ws.gf_ref, ws.fname)
ws.Compare(ws.gf_regridded, ws.gf_ref, ws.maxabsdiff)
#
# Maximum pressure in new grid is higher than in the old grid
#
ws.StringSet(ws.fname, "gfREFERENCE_z_0p_bottom.xml")
ws.ReadXML(ws.p_grid, "p_grid_bottom.xml")
ws.ReadXML(ws.z_field_raw, "z_bottom.xml")
ws.FieldFromGriddedField(
    ws.z_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.z_field_raw
)
ws.GriddedFieldZToPRegrid(
    ws.gf_regridded, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.z_field, ws.gf_z, 1, 1
)
# WriteXML("ascii", gf_regridded, fname)
ws.ReadXML(ws.gf_ref, ws.fname)
ws.Compare(ws.gf_regridded, ws.gf_ref, ws.maxabsdiff)
#
# Minimum pressure in new grid is lower than in the old grid
#
ws.StringSet(ws.fname, "gfREFERENCE_z_0p_top.xml")
ws.ReadXML(ws.p_grid, "p_grid_top.xml")
ws.ReadXML(ws.z_field_raw, "z_top.xml")
ws.FieldFromGriddedField(
    ws.z_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.z_field_raw
)
ws.GriddedFieldZToPRegrid(
    ws.gf_regridded, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.z_field, ws.gf_z, 1, 1
)
# WriteXML("ascii", gf_regridded, fname)
ws.ReadXML(ws.gf_ref, ws.fname)
ws.Compare(ws.gf_regridded, ws.gf_ref, ws.maxabsdiff)
