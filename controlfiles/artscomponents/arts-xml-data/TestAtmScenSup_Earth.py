#
# Testing functionality (meeting format requirements, etc.) of supplemental
#  atmospheric scenario data.
#
# General test setup: reading in raw data (including a basic atmosphere),
#  expanding (where necessary), regridding (incl. extracting), executing
#  standard pre-RT calc internal test method atmfields_checkedCalc.
#
#
# This case is for Earth and specifically tests
#
#  - electron densities (given: 4seasons x 8dayttimes x 2solaractivity global 3D
#    cases)
#      - all 64 cases in global 3D with p_grid from altitude grid taken from
#        standard FASCOD, i.e. extending up to 95km only (CASEs A) and p_grid
#        from altitude grid taken from higher altitude extensions (up to 2000km)
#        to FASCOD (CASEs B)
#  - magnetic field (given: 4decadal 3D cases, each with with 3 components)
#      - global 3D case with p_grid from standard and expanded FASCOD altitude
#        grids (CASEs A and B)
#
# NOTE: beside applying them here as "background" atmosphere, the FASCOD
#  atmospheric field data is not tested. This dataset has been part of ARTS
#  xml-data for long time, applied in calculations, hence further testing is not
#  necessary.
#
# Jana Mendrok 2013-02-26

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
ws.Tensor3Create("edensity_field")
ws.GriddedField3Create("gf3tmp")
ws.StringCreate("caseext")
ws.StringCreate("casefull")
ws.StringCreate("basename")
ws.StringCreate("atmcase")
ws.IndexCreate("ncases")
ws.IndexCreate("interp_order")
ws.IndexSet(ws.interp_order, 1)
# Array with base case names
ws.ArrayOfStringCreate("basecasearray")
ws.ArrayOfStringSet(
    ws.basecasearray,
    ["planets/Earth/Fascod/", "planets/Earth/IRI/IRI", "planets/Earth/IGRF/IGRF11"],
)
#   "/storage4/home/mendrok/projects/MicrowavePropagationToolbox/WorkPackageData/WP2000/Earth/MagField/ToolBoxdata/IGRF11"] )
# StringSet( atmcase, "midlatitude-summer" )
ws.StringSet(ws.atmcase, "tropical")
# Arrays with Ne (sub)case names
ws.ArrayOfStringCreate("seasoncasearray")
ws.ArrayOfStringSet(ws.seasoncasearray, ["_spring", "_summer", "_fall", "_winter"])
ws.ArrayOfStringCreate("timecasearray")
ws.ArrayOfStringSet(
    ws.timecasearray,
    [
        "_00UTC.Ne",
        "_03UTC.Ne",
        "_06UTC.Ne",
        "_09UTC.Ne",
        "_12UTC.Ne",
        "_15UTC.Ne",
        "_18UTC.Ne",
        "_21UTC.Ne",
    ],
)
ws.ArrayOfStringCreate("solarcasearray")
ws.ArrayOfStringSet(ws.solarcasearray, ["_solmax", "_solmin"])
# Arrays with B (sub)case names
ws.ArrayOfStringCreate("yearcasearray")
ws.ArrayOfStringSet(ws.yearcasearray, ["_1980", "_1990", "_2000", "_2010"])
# , "_2012"] )
ws.ArrayOfStringCreate("resolcasearray")
ws.ArrayOfStringSet(ws.resolcasearray, ["_200km-5deg-5deg"])
# the 3D geo grid to test
ws.VectorCreate("lat_grid3D")
ws.VectorCreate("lon_grid3D")
ws.VectorLinSpace(ws.lat_grid3D, -90.0, 90.0, 18.0)
ws.VectorLinSpace(ws.lon_grid3D, -20.0, 350.0, 18.0)
ws.VectorLinSpace(ws.lon_grid3D, 0.0, 360.0, 18.0)
#####
# first electron densities, 3D
#####
###
# this is the real stuff, partI: electron density
# We have to repeat that, as we have multiple Ne fields given, but each atm case
# can have only one (in contrast to vmr profiles, where we could test all
# profiles in a single run). We use a bunch of forloops for this, which we
# define here and execute later on.
ws.AgendaCreate("forloop_agenda_time")


@arts_agenda
def forloop_agenda_time(ws):
    # construct case name
    ws.Extract(ws.casefull, ws.timecasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Print(ws.atmcase, 0)
    # readin in raw data
    ws.ReadXML(ws.gf3tmp, ws.atmcase)
    # this is 3D data. we need to regrid the raw data to the calculation grid(s).
    # for supplemental atm data (Ne, B, winds) this requires manual regridding (in
    # contrast to basic atm data, where this is handled by AtmFieldsCalc(Expand1D).
    ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
    ws.GriddedFieldLatLonRegrid(
        ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
    )
    ws.FieldFromGriddedField(
        ws.edensity_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
    )
    ws.atmfields_checkedCalc()
    ws.atmgeom_checkedCalc()
    # WriteXML( "ascii", p_grid )
    # WriteXML( "ascii", z_field )
    # WriteXML( "ascii", t_field )
    # WriteXMLIndexed( "ascii", forloop_index, edensity_field )


ws.forloop_agenda_time = forloop_agenda_time

ws.AgendaCreate("forloop_agenda_season")


@arts_agenda
def forloop_agenda_season(ws):
    # construct case name
    ws.Extract(ws.casefull, ws.seasoncasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    # Print( atmcase, 0 )
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_time)
    ws.nelemGet(ws.ncases, ws.timecasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_season = forloop_agenda_season

ws.AgendaCreate("forloop_agenda_solar")


@arts_agenda
def forloop_agenda_solar(ws):
    # construct case name
    ws.Copy(ws.atmcase, ws.basename)
    ws.Extract(ws.casefull, ws.solarcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    # Print( atmcase, 0 )
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_season)
    ws.nelemGet(ws.ncases, ws.seasoncasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_solar = forloop_agenda_solar

# now we actually execute things...
# 3-dimensional atmosphere
ws.AtmosphereSet3D()
ws.Copy(ws.lat_grid, ws.lat_grid3D)
ws.Copy(ws.lon_grid, ws.lon_grid3D)
# set atmospheric scenario
ws.Extract(ws.basename, ws.basecasearray, 0)
ws.Append(ws.basename, ws.atmcase)
ws.StringSet(ws.caseext, "/")
ws.Append(ws.basename, ws.caseext)
ws.Append(ws.basename, ws.atmcase)
# we manually select a minumim set of basic atm data (main atm constituents)
ws.abs_speciesSet(species=["O2", "N2"])
ws.AtmRawRead(basename=ws.basename)
ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
ws.AtmFieldsCalcExpand1D(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 0)
# now doing the Ne cases (reading, 1D regridding, checking). we use several
# forloops to loop over the different aspects of the data
ws.Extract(ws.basename, ws.basecasearray, 1)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_solar)
ws.nelemGet(ws.ncases, ws.solarcasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
# now we want to use a larger p_grid (for using the full vertical extend of the
# Ne data). that means, we need to get larger extending z_field and t_field. we
# can derive a z_field from p_grid and t_field, but we need a larger t_field. as
# manipulating fields within ARTS controlfiles is a hassle, we made extended raw
# t/z fields that as far as possible agree with the assumptions made when
# creating the Ne and B fields (IRI and IGRF data). we replace the standard t/z
# raw field by these.
# going back to clearsky atmosphere and derive expanded t/z
ws.Extract(ws.basename, ws.basecasearray, 0)
ws.Append(ws.basename, ws.atmcase)
ws.StringSet(ws.caseext, "/")
ws.Append(ws.basename, ws.caseext)
ws.Append(ws.basename, ws.atmcase)
ws.Copy(ws.casefull, ws.basename)
ws.StringSet(ws.caseext, ".expanded.t")
ws.Append(ws.casefull, ws.caseext)
ws.Print(ws.casefull, 0)
ws.ReadXML(ws.t_field_raw, ws.casefull)
ws.Copy(ws.casefull, ws.basename)
ws.StringSet(ws.caseext, ".expanded.z")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.z_field_raw, ws.casefull)
ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
# make a higher resolution p_grid spread over input z_field (i.e.,with same
# upper and lower boundaries, but more grid points)
ws.IndexCreate("itmp")
ws.NumericCreate("ntmp1")
ws.NumericCreate("ntmp2")
ws.Extract(ws.ntmp1, ws.p_grid, 0)
ws.nelemGet(ws.itmp, ws.p_grid)
ws.IndexStepDown(ws.itmp, ws.itmp)
ws.Extract(ws.ntmp2, ws.p_grid, ws.itmp)
ws.VectorNLogSpace(ws.p_grid, 201, ws.ntmp1, ws.ntmp2)
# WriteXML( "ascii", p_grid )
ws.AtmFieldsCalcExpand1D(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 0)
# now doing the Ne cases (reading, 1D regridding, checking). we use several
# forloops to loop over the different aspects of the data
ws.Extract(ws.basename, ws.basecasearray, 1)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_solar)
ws.nelemGet(ws.ncases, ws.solarcasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
#####
# second: magnetic field, 3D
#####
ws.AgendaCreate("forloop_agenda_Bres")


@arts_agenda
def forloop_agenda_Bres(ws):
    # construct case name
    ws.Extract(ws.casefull, ws.resolcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Print(ws.atmcase, 0)
    # this is 3D data. we need to regrid the raw data to the calculation grid(s).
    # for supplemental atm data (Ne, B, winds) this requires manual reading in and
    # regridding (in contrast to basic atm data, where this is handled by
    # AtmFieldsCalc(Expand1D).
    # now reading and regridding of each of the 3 B-field components separately
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".mag_u")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
    ws.GriddedFieldLatLonRegrid(
        ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
    )
    ws.FieldFromGriddedField(
        ws.mag_u_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
    )
    # now reading and regridding of each of the 3 B-field components separately
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".mag_v")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
    ws.GriddedFieldLatLonRegrid(
        ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
    )
    ws.FieldFromGriddedField(
        ws.mag_v_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
    )
    # now reading and regridding of each of the 3 B-field components separately
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".mag_w")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
    ws.GriddedFieldLatLonRegrid(
        ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
    )
    ws.FieldFromGriddedField(
        ws.mag_w_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
    )
    ws.atmfields_checkedCalc()
    ws.atmgeom_checkedCalc()
    # WriteXML( "ascii", p_grid )
    # WriteXML( "ascii", z_field )
    # WriteXML( "ascii", t_field )
    # WriteXMLIndexed( "ascii", forloop_index, mag_u_field )
    # WriteXMLIndexed( "ascii", forloop_index, mag_v_field )
    # WriteXMLIndexed( "ascii", forloop_index, mag_w_field )


ws.forloop_agenda_Bres = forloop_agenda_Bres

ws.AgendaCreate("forloop_agenda_Byears")


@arts_agenda
def forloop_agenda_Byears(ws):
    # construct case name
    ws.Copy(ws.atmcase, ws.basename)
    ws.Extract(ws.casefull, ws.yearcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Print(ws.atmcase, 0)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_Bres)
    ws.nelemGet(ws.ncases, ws.resolcasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_Byears = forloop_agenda_Byears

# now we actually execute things...
# 3-dimensional atmosphere
ws.AtmosphereSet3D()
ws.Copy(ws.lat_grid, ws.lat_grid3D)
ws.Copy(ws.lon_grid, ws.lon_grid3D)
# set atmospheric scenario
ws.Extract(ws.basename, ws.basecasearray, 0)
ws.Append(ws.basename, ws.atmcase)
ws.StringSet(ws.caseext, "/")
ws.Append(ws.basename, ws.caseext)
ws.Append(ws.basename, ws.atmcase)
# we manually select a minumim set of basic atm data (main atm constituents)
ws.Delete(ws.vmr_field_raw)
ws.abs_speciesSet(species=["O2", "N2"])
ws.AtmRawRead(basename=ws.basename)
ws.Copy(ws.casefull, ws.basename)
ws.StringSet(ws.caseext, ".expanded.t")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.t_field_raw, ws.casefull)
ws.Copy(ws.casefull, ws.basename)
ws.StringSet(ws.caseext, ".expanded.z")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.z_field_raw, ws.casefull)
ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
# make a higher resolution p_grid spread over input z_field (i.e.,with same
# upper and lower boundaries, but more grid points)
ws.Extract(ws.ntmp1, ws.p_grid, 0)
ws.nelemGet(ws.itmp, ws.p_grid)
ws.IndexStepDown(ws.itmp, ws.itmp)
ws.Extract(ws.ntmp2, ws.p_grid, ws.itmp)
ws.VectorNLogSpace(ws.p_grid, 201, ws.ntmp1, ws.ntmp2)
# WriteXML( "ascii", p_grid )
ws.AtmFieldsCalcExpand1D(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
# now doing the B cases (reading, 1D regridding, checking). we use several
# forloops to loop over the different aspects of the data
ws.Extract(ws.basename, ws.basecasearray, 2)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_Byears)
ws.nelemGet(ws.ncases, ws.yearcasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
