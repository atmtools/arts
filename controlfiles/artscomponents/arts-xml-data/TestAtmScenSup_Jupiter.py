#
# Testing functionality (meeting format requirements, etc.) of supplemental
#  atmospheric scenario data.
#
# General test setup: reading in raw data (including a basic atmosphere),
#  expanding (where necessary), regridding (incl. extracting), executing
#  standard pre-RT calc internal test method atmfields_checkedCalc.
#
#
# This case is for Jupiter and specifically tests
#
#  - electron densities (given: 2x3 1D cases)
#      - all cases from Jupiter.mean and Jupiter.oval (each with 3 density
#        levels) in 1D with p_grid from altitude grid (CASE A) and user defined
#        p_grid (CASE B)
#      - a single case expanded to 3D (assuming the other cases behave in the
#        same way) with p_grid from altitude grid
#  - magnetic field (given: single 3D case with 3 components)
#      - extracting 1D case from 3D with p_grid from altitude grid
#      - global 3D case with p_grid from altitude grid (also tested with B's own
#        z_field)
#  - wind (given: 2 1D cases with 1 component)
#      - both expanded to global 3D case with p_grid from altitude grid
#
# Jana Mendrok 2013-02-26

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_jupiter.arts")
# 1-dimensional atmosphere
ws.AtmosphereSet1D()
# We have to set the absorption lookup table interpolation order,
# since we are using wind fields. atmfields_checkedCalc will otherwise throw an error.
ws.IndexSet(ws.abs_f_interp_order, 1)
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
    [
        "planets/Jupiter/MPS/Jupiter.mean/Jupiter.mean",
        "planets/Jupiter/MPS/Jupiter.oval/Jupiter.oval",
        "planets/Jupiter/Khurana/Khurana",
    ],
)
# Array with Ne case names
ws.ArrayOfStringCreate("necasearray")
ws.ArrayOfStringSet(ws.necasearray, [".Ne.high", ".Ne.low", ".Ne.med"])
###
# this is the real stuff, partI: electron density
# We have to repeat that, as we have multiple Ne fields given, but each atm case
# can have only one (in contrast to vmr profiles, where we could test all
# profiles in a single run). We use a forloop for this, which we define here and
# execute later on.
ws.AgendaCreate("forloop_agenda_Ne")


@arts_agenda
def forloop_agenda_Ne(ws):
    # construct case name
    ws.Copy(ws.atmcase, ws.basename)
    ws.Extract(ws.casefull, ws.necasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Print(ws.atmcase, 0)
    # readin in raw data
    ws.ReadXML(ws.gf3tmp, ws.atmcase)
    # this is 1D data and we're doing 1D. but we need to regrid the raw data to
    # the calculation grid(s). for supplemental atm data (Ne, B, winds) this
    # requires manual regridding (in contrast to basic atm data, where this is
    # handled by AtmFieldsCalc.
    # so, first: regrid to p_grid (as we are in 1D, we don't need latlon regridding)
    ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
    # eventually, extract the data Tensor from the regridded GriddedField
    ws.FieldFromGriddedField(
        ws.edensity_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
    )
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()
    # WriteXML( "ascii", p_grid )
    # WriteXML( "ascii", z_field )
    # WriteXML( "ascii", t_field )
    # WriteXMLIndexed( "ascii", forloop_index, edensity_field )


ws.forloop_agenda_Ne = forloop_agenda_Ne

#####
# Ne: CASE A
#####
# set atmospheric scenario
ws.Extract(ws.basename, ws.basecasearray, 0)
# we manually select a minumim set of basic atm data (main atm constituents)
ws.abs_speciesSet(species=["H2", "He"])
# get atm scenario raw data
ws.AtmRawRead(basename=ws.basename)
#####
# A-1) p_grid initialized from given altitude grid
#####
# we need to include negative altitudes!
ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
ws.AtmFieldsCalc(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 0)
# now looping over the Ne cases (reading, 1D regridding, checking)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_Ne)
ws.nelemGet(ws.ncases, ws.necasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
#####
# A-2) p_grid set to a user defined grid (surely requries interpolation to calc-grid(s))
#####
ws.VectorNLogSpace(ws.p_grid, 401, 920000.0, 0.001)
ws.AtmFieldsCalc(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 0)
# repeating forloop with different p_grid
ws.nelemGet(ws.ncases, ws.necasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
#####
# Ne: CASE B
#####
# Now Ne from the second case folder. For the rest, we use the profiles (t, z,
# vmr) from above.
# set atmospheric scenario
ws.Extract(ws.basename, ws.basecasearray, 1)
# Copy( casefull, basename )
# StringSet( caseext, ".t" )
# Append( casefull, caseext )
# ReadXML( t_field_raw, casefull )
# Copy( casefull, basename )
# StringSet( caseext, ".z" )
# Append( casefull, caseext )
# ReadXML( z_field_raw, casefull )
#####
# B-1) p_grid initialized from given altitude grid
#####
# we need to include negative altitudes!
ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
ws.AtmFieldsCalc(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 0)
# repeating forloop with CASE B Ne scenario
ws.nelemGet(ws.ncases, ws.necasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
#####
# Magfield
#####
# We continue with magnetic field (B). This is given in 3D, but first we try a 1D
# case, i.e. we will extract data for a given geographical location.
# We have only one set of B-field data given, hence we don't need a forloop here.
ws.VectorSet(ws.lat_true, array([-12.0]))
ws.VectorSet(ws.lon_true, array([10.0]))
# construct atmcase
ws.Extract(ws.basename, ws.basecasearray, 2)
# now reading and regridding of each of the 3 B-field components separately
ws.Copy(ws.atmcase, ws.basename)
ws.StringSet(ws.caseext, ".B_u")
ws.Append(ws.atmcase, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.atmcase)
# first: regrid to p_grid
ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
# second: regrid to lat/lon_grid
ws.GriddedFieldLatLonRegrid(
    ws.gf3tmp, ws.lat_true, ws.lon_true, ws.gf3tmp, ws.interp_order
)
# last, extract the data Tensor from the regridded GriddedField
ws.FieldFromGriddedField(ws.mag_u_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp)
ws.Copy(ws.atmcase, ws.basename)
ws.StringSet(ws.caseext, ".B_v")
ws.Append(ws.atmcase, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.atmcase)
# first: regrid to p_grid
ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
# second: regrid to lat/lon_grid
ws.GriddedFieldLatLonRegrid(
    ws.gf3tmp, ws.lat_true, ws.lon_true, ws.gf3tmp, ws.interp_order
)
# last, extract the data Tensor from the regridded GriddedField
ws.FieldFromGriddedField(ws.mag_v_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp)
ws.Copy(ws.atmcase, ws.basename)
ws.StringSet(ws.caseext, ".B_w")
ws.Append(ws.atmcase, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.atmcase)
# first: regrid to p_grid
ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
# second: regrid to lat/lon_grid
ws.GriddedFieldLatLonRegrid(
    ws.gf3tmp, ws.lat_true, ws.lon_true, ws.gf3tmp, ws.interp_order
)
# last, extract the data Tensor from the regridded GriddedField
ws.FieldFromGriddedField(ws.mag_w_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp)
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
# WriteXML( "ascii", mag_u_field )
# WriteXML( "ascii", mag_v_field )
# WriteXML( "ascii", mag_w_field )
#####
# and now in 3D
#####
# grid settings
ws.AtmosphereSet3D()
ws.VectorLinSpace(ws.lat_grid, -90.0, 90.0, 18.0)
ws.VectorLinSpace(ws.lon_grid, -20.0, 340.0, 18.0)
# blowing up basic atmosphere (and surface)
ws.AtmFieldsExpand1D()
ws.Extract(ws.z_surface, ws.z_field, 0)
# rereading (we need the GriddedField, which we didn't keep), blowing up, and
# regridding Ne. but we take only one of them, and assume the others to behave
# in the same way.
ws.Extract(ws.atmcase, ws.basecasearray, 1)
ws.Extract(ws.casefull, ws.necasearray, 2)
ws.Append(ws.atmcase, ws.casefull)
ws.Print(ws.atmcase, 0)
ws.ReadXML(ws.gf3tmp, ws.atmcase)
ws.GriddedFieldLatLonExpand(ws.gf3tmp, ws.gf3tmp)
ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
ws.GriddedFieldLatLonRegrid(
    ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
)
ws.FieldFromGriddedField(
    ws.edensity_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
)
# rereading (since we didn't keep the 3D raw data) and regridding B-field
# construct atmcase
ws.Extract(ws.basename, ws.basecasearray, 2)
# now reading and regridding of each of the 3 B-field components separately
ws.Copy(ws.atmcase, ws.basename)
ws.StringSet(ws.caseext, ".B_u")
ws.Append(ws.atmcase, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.atmcase)
# first: regrid to p_grid
ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
# second: regrid to lat/lon_grid
ws.GriddedFieldLatLonRegrid(
    ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
)
# last, extract the data Tensor from the regridded GriddedField
ws.FieldFromGriddedField(ws.mag_u_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp)
ws.Copy(ws.atmcase, ws.basename)
ws.StringSet(ws.caseext, ".B_v")
ws.Append(ws.atmcase, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.atmcase)
# first: regrid to p_grid
ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
# second: regrid to lat/lon_grid
ws.GriddedFieldLatLonRegrid(
    ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
)
# last, extract the data Tensor from the regridded GriddedField
ws.FieldFromGriddedField(ws.mag_v_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp)
ws.Copy(ws.atmcase, ws.basename)
ws.StringSet(ws.caseext, ".B_w")
ws.Append(ws.atmcase, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.atmcase)
# first: regrid to p_grid
ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
# second: regrid to lat/lon_grid
ws.GriddedFieldLatLonRegrid(
    ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
)
# last, extract the data Tensor from the regridded GriddedField
ws.FieldFromGriddedField(ws.mag_w_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp)
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
#####
# Winds
#####
# reading and regridding of the wind-field components. however, here we only have
# data for u-component (zonal), but two different version of them.
ws.Extract(ws.basename, ws.basecasearray, 0)
###
# version 1
ws.Copy(ws.atmcase, ws.basename)
ws.StringSet(ws.caseext, ".wind_u")
ws.Append(ws.atmcase, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.atmcase)
# first: regrid to p_grid
ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
# second: make 3D fields from 1D, then regrid to lat/lon_grid
ws.GriddedFieldLatLonExpand(ws.gf3tmp, ws.gf3tmp)
ws.GriddedFieldLatLonRegrid(
    ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
)
# last, extract the data Tensor from the regridded GriddedField
ws.FieldFromGriddedField(
    ws.wind_u_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
)
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
###
# version 2
ws.Copy(ws.atmcase, ws.basename)
ws.StringSet(ws.caseext, ".wind_u.into-thermal")
ws.Append(ws.atmcase, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.atmcase)
# first: regrid to p_grid
ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
# second: make 3D fields from 1D, then regrid to lat/lon_grid
ws.GriddedFieldLatLonExpand(ws.gf3tmp, ws.gf3tmp)
ws.GriddedFieldLatLonRegrid(
    ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
)
# last, extract the data Tensor from the regridded GriddedField
ws.FieldFromGriddedField(
    ws.wind_u_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
)
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
# WriteXML( "ascii", wind_u_field )
# WriteXML( "ascii", wind_v_field )
# WriteXML( "ascii", wind_w_field )
