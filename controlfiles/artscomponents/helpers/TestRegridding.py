#
# Testing functionality of Regridding routines GriddedFieldLatLonRegrid and
#  GriddedFieldPRegrid
#
# Jana Mendrok 2013-02-26

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.StringCreate("atmcase")
ws.StringSet(ws.atmcase, "planets/Earth/Fascod/tropical/tropical")
ws.StringCreate("surfpath")
ws.StringSet(ws.surfpath, "planets/Earth/ECMWF/ERA40/")
ws.StringCreate("zsurfname")
ws.StringSet(ws.zsurfname, "SurfaceAltitude_ERA40_1.0Degree")
ws.StringCreate("surfmaskname")
ws.StringSet(ws.surfmaskname, "LandSeaMask_ERA40_1.0Degree")
ws.StringCreate("Bname")
ws.StringSet(ws.Bname, "planets/Earth/IGRF/IGRF11_2010_200km-5deg-5deg")
ws.StringCreate("windname")
ws.StringSet(
    ws.windname,
    "planets/Mars/MPS/Mars.Ls0.day.dust-high/Mars.Ls0.day.dust-high.sol-avg/Mars.Ls0.day.dust-high.sol-avg",
)
ws.StringCreate("casefull")
ws.StringCreate("caseext")
ws.StringCreate("slat")
ws.StringSet(ws.slat, "lat:")
ws.StringCreate("slon")
ws.StringSet(ws.slon, "lon:")
ws.VectorNLogSpace(ws.p_grid, 11, 100000.0, 0.1)
ws.VectorCreate("lat_regrid")
ws.VectorLinSpace(ws.lat_regrid, 45.0, -20.0, 10.0)
ws.VectorCreate("lon_regrid")
ws.VectorLinSpace(ws.lon_regrid, -3.5, 12.5, 0.25)
ws.NumericCreate("lattrue")
ws.NumericCreate("lontrue")
ws.IndexCreate("ncases")
ws.IndexCreate("interp_order")
ws.IndexSet(ws.interp_order, 1)
ws.GriddedField2Create("z_surface_raw")
ws.GriddedField2Create("lsmask_raw")
ws.MatrixCreate("lsmask")
ws.GriddedField3Create("B_field_raw")
ws.GriddedField3Create("wind_field_raw")
ws.MatrixCreate("mtmp")
ws.Tensor3Create("t3tmp")
ws.abs_speciesSet(species=["CH4", "H2O", "O2", "N2"])
# get atm scenario raw data
ws.AtmRawRead(basename=ws.atmcase)
# for testing and demonstration purposes, we manually expand 1D atmosphere first
# to 3D then regrid back to 1D
ws.GriddedFieldLatLonExpand(ws.t_field_raw, ws.t_field_raw)
ws.GriddedFieldLatLonExpand(ws.z_field_raw, ws.z_field_raw)
ws.GriddedFieldLatLonExpand(ws.vmr_field_raw, ws.vmr_field_raw)
# reading the surface altitude field
ws.Copy(ws.casefull, ws.surfpath)
ws.Append(ws.casefull, ws.zsurfname)
ws.ReadXML(ws.z_surface_raw, ws.casefull)
# reading surface mask
ws.Copy(ws.casefull, ws.surfpath)
ws.Append(ws.casefull, ws.surfmaskname)
ws.ReadXML(ws.lsmask_raw, ws.casefull)
# reading B-field w-component
ws.Copy(ws.casefull, ws.Bname)
ws.StringSet(ws.caseext, ".mag_w.xml")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.B_field_raw, ws.casefull)
# reading wind-field w-component
ws.Copy(ws.casefull, ws.windname)
ws.StringSet(ws.caseext, ".wind_w.xml")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.wind_field_raw, ws.casefull)
ws.AgendaCreate("forloop_agenda_lon")


@arts_agenda
def forloop_agenda_lon(ws):
    ws.Extract(ws.lontrue, ws.lon_regrid, ws.forloop_index)
    ws.VectorSetConstant(ws.lon_true, 1, ws.lontrue)
    # Print( slon, 0 )
    # Print( lon_true, 0 )
    ws.GriddedFieldLatLonRegrid(
        ws.t_field_raw, ws.lat_true, ws.lon_true, ws.t_field_raw
    )
    ws.GriddedFieldPRegrid(
        ws.t_field_raw, ws.p_grid, ws.t_field_raw, ws.interp_order, 0
    )
    ws.FieldFromGriddedField(
        ws.t_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.t_field_raw
    )
    ws.Extract(ws.mtmp, ws.t_field, 5)
    # Print( mtmp, 0 )
    ws.GriddedFieldLatLonRegrid(
        ws.z_field_raw, ws.lat_true, ws.lon_true, ws.z_field_raw
    )
    ws.GriddedFieldPRegrid(
        ws.z_field_raw, ws.p_grid, ws.z_field_raw, ws.interp_order, 0
    )
    ws.FieldFromGriddedField(
        ws.z_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.z_field_raw
    )
    ws.Extract(ws.mtmp, ws.z_field, 5)
    # Print( mtmp, 0 )
    ws.GriddedFieldLatLonRegrid(
        ws.vmr_field_raw, ws.lat_true, ws.lon_true, ws.vmr_field_raw
    )
    ws.GriddedFieldPRegrid(
        ws.vmr_field_raw, ws.p_grid, ws.vmr_field_raw, ws.interp_order, 1
    )
    ws.FieldFromGriddedField(
        ws.vmr_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.vmr_field_raw
    )
    ws.Extract(ws.t3tmp, ws.vmr_field, 0)
    ws.Extract(ws.mtmp, ws.t3tmp, 5)
    # Print( mtmp, 0 )
    ws.GriddedFieldLatLonRegrid(
        ws.z_surface_raw, ws.lat_true, ws.lon_true, ws.z_surface_raw
    )
    ws.FieldFromGriddedField(
        ws.z_surface, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.z_surface_raw
    )
    # Print( z_surface, 0 )
    ws.GriddedFieldLatLonRegrid(ws.lsmask_raw, ws.lat_true, ws.lon_true, ws.lsmask_raw)
    ws.FieldFromGriddedField(
        ws.lsmask, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.lsmask_raw
    )
    # Print( lsmask, 0 )
    ws.GriddedFieldPRegrid(
        ws.B_field_raw, ws.p_grid, ws.B_field_raw, ws.interp_order, 1
    )
    ws.GriddedFieldLatLonRegrid(
        ws.B_field_raw, ws.lat_true, ws.lon_true, ws.B_field_raw, ws.interp_order
    )
    ws.FieldFromGriddedField(
        ws.mag_w_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.B_field_raw
    )
    ws.Extract(ws.mtmp, ws.mag_w_field, 5)
    # Print( mtmp, 0 )
    # we can also do the expanding here...
    ws.GriddedFieldPRegrid(
        ws.wind_field_raw, ws.p_grid, ws.wind_field_raw, ws.interp_order, 1
    )
    ws.GriddedFieldLatLonExpand(ws.wind_field_raw, ws.wind_field_raw)
    ws.GriddedFieldLatLonRegrid(
        ws.wind_field_raw, ws.lat_true, ws.lon_true, ws.wind_field_raw, ws.interp_order
    )
    ws.FieldFromGriddedField(
        ws.wind_w_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.wind_field_raw
    )
    ws.Extract(ws.mtmp, ws.wind_w_field, 5)
    # Print( mtmp, 0 )


ws.forloop_agenda_lon = forloop_agenda_lon

ws.AgendaCreate("forloop_agenda_lat")


@arts_agenda
def forloop_agenda_lat(ws):
    ws.Extract(ws.lattrue, ws.lat_regrid, ws.forloop_index)
    ws.VectorSetConstant(ws.lat_true, 1, ws.lattrue)
    # Print( slat, 0 )
    # Print( lat_true, 0 )
    ws.nelemGet(ws.ncases, ws.lon_regrid)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_lon)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_lat = forloop_agenda_lat

ws.AtmosphereSet1D()
ws.nelemGet(ws.ncases, ws.lat_regrid)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_lat)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
