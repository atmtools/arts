#
# Testing functionality (meeting format requirements, etc.) of supplemental
#  atmospheric scenario data.
#
# General test setup: reading in raw data (including a basic atmosphere),
#  expanding (where necessary), regridding (incl. extracting), executing
#  standard pre-RT calc internal test method atmfields_checkedCalc.
#
#
# This case is for Mars and specifically tests
#
#  - for the 72 Mars scenarios: 4seasons x 2daytimes x 3dustloads x 3solaractivities
#  - electron densities (given: 1D cases, 1 per night scenario, 4 per day scenario)
#      - all cases in 1D with p_grid from altitude grid
#      - a single case per scenario expanded to 3D (assuming the other cases
#        behave in the same way) with p_grid from altitude grid
#  - wind (given: one 1D case per scenario with all 3 components)
#      - all expanded to global 3D case with p_grid from altitude grid
#
# Jana Mendrok 2013-02-26

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_mars.arts")
# 1-dimensional atmosphere
ws.AtmosphereSet1D()
# We have to set the absorption lookup table interpolation order,
# since we are using wind fields. atmfields_checkedCalc will otherwise throw an error.
ws.IndexSet(ws.abs_f_interp_order, 1)
ws.Tensor3Create("edensity_field")
ws.GriddedField3Create("gf3tmp")
ws.StringCreate("caseext")
ws.StringCreate("casefull")
ws.StringCreate("atmcase")
ws.IndexCreate("ncases")
ws.IndexCreate("interp_order")
ws.IndexSet(ws.interp_order, 1)
# set basic case folder
ws.StringCreate("basename")
ws.StringSet(ws.basename, "planets/Mars/MPS/")
# Arrays with (sub)case names
ws.ArrayOfStringCreate("seasoncasearray")
ws.ArrayOfStringSet(
    ws.seasoncasearray, ["Mars.Ls0", "Mars.Ls90", "Mars.Ls180", "Mars.Ls270"]
)
ws.ArrayOfStringCreate("timecasearray")
ws.ArrayOfStringSet(ws.timecasearray, [".day", ".night"])
ws.ArrayOfStringCreate("dustcasearray")
ws.ArrayOfStringSet(ws.dustcasearray, [".dust-high", ".dust-low", ".dust-medium"])
ws.ArrayOfStringCreate("solarcasearray")
ws.ArrayOfStringSet(ws.solarcasearray, [".sol-avg", ".sol-max", ".sol-min"])
ws.ArrayOfStringCreate("necasearray")
ws.ArrayOfStringCreate("nightnearray")
ws.ArrayOfStringSet(ws.nightnearray, [".SZA.120-180.Ne"])
ws.ArrayOfStringCreate("daynearray")
ws.ArrayOfStringSet(
    ws.daynearray, [".SZA.0-30.Ne", ".SZA.30-50.Ne", ".SZA.50-70.Ne", ".SZA.70-90.Ne"]
)
# the 3D geo grid to test
ws.VectorCreate("lat_grid3D")
ws.VectorCreate("lon_grid3D")
ws.VectorLinSpace(ws.lat_grid3D, -90.0, 90.0, 18.0)
ws.VectorLinSpace(ws.lon_grid3D, -20.0, 350.0, 18.0)
# we go with several nested foorloop through the different cases.
#  All those cases have identical abs species to process.
#  Order of agenda definitions has to be inverse from their execution (as we can
#  only copy an agenda AFTER we have defined it).
ws.AgendaCreate("forloop_agenda_Ne")


@arts_agenda
def forloop_agenda_Ne(ws):
    # construct the full case name
    ws.Copy(ws.casefull, ws.basename)
    ws.Extract(ws.caseext, ws.necasearray, ws.forloop_index)
    ws.Append(ws.casefull, ws.caseext)
    # Print( casefull, 0 )
    # readin in raw data
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    # this is 1D data and we're doing 1D. but we need to regrid the raw data to
    # the calculation grid(s). for supplemental atm data (Ne, B, winds) this
    # requires manual regridding (in contrast to basic atm data, where this is
    # handled by AtmFieldsCalc.
    # so, first: regrid to p_grid (as we are in 1D, we don't need latlon regridding)
    ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
    # eventually: extract the data Tensor from the regridded GriddedField
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

ws.AgendaCreate("forloop_agenda_solar")


@arts_agenda
def forloop_agenda_solar(ws):
    # construct atmcase name IV (Mars.LsXX.YY.dust-ZZ.sol-WW)
    ws.Extract(ws.casefull, ws.solarcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Append(ws.basename, ws.atmcase)
    ws.StringSet(ws.caseext, "/")
    ws.Append(ws.basename, ws.caseext)
    ws.Append(ws.basename, ws.atmcase)
    # Print( basename, 0 )
    # we manually select a minumim set of basic atm data (main atm constituents)
    ws.abs_speciesSet(species=["CO2"])
    ws.AtmRawRead(basename=ws.basename)
    ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
    ws.AtmFieldsCalc()
    ws.Extract(ws.z_surface, ws.z_field, 0)
    # now get the Ne data. with several cases per scenario, we use another forloop
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_Ne)
    ws.nelemGet(ws.ncases, ws.necasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
    # now changing to 3D for the winds (but do one Ne case per scenario as well)
    ws.AtmosphereSet3D()
    ws.Copy(ws.lat_grid, ws.lat_grid3D)
    ws.Copy(ws.lon_grid, ws.lon_grid3D)
    # blowing up basic atmosphere (and surface)
    ws.AtmFieldsExpand1D()
    ws.Extract(ws.z_surface, ws.z_field, 0)
    # getting and preprocessing Ne data in 3D (single case per scenario)
    ws.Copy(ws.casefull, ws.basename)
    ws.Extract(ws.caseext, ws.necasearray, 0)
    ws.Append(ws.casefull, ws.caseext)
    # Print( casefull, 0 )
    # reading in raw data
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.GriddedFieldLatLonExpand(ws.gf3tmp, ws.gf3tmp)
    ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
    ws.GriddedFieldLatLonRegrid(
        ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
    )
    ws.FieldFromGriddedField(
        ws.edensity_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
    )
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()
    # and eventually the winds - here we the two horizontal components (u&v), and
    # the u-component (zonal) in 3 different cases.
    # first the w-component
    ws.Copy(ws.casefull, ws.basename)
    ws.StringSet(ws.caseext, ".wind_w")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    # Print( casefull, 0 )
    # first: regrid to p_grid
    ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
    # second: make 3D fields from 1D, then regrid to lat/lon_grid
    ws.GriddedFieldLatLonExpand(ws.gf3tmp, ws.gf3tmp)
    ws.GriddedFieldLatLonRegrid(
        ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
    )
    # last, extract the data Tensor from the regridded GriddedField
    ws.FieldFromGriddedField(
        ws.wind_w_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
    )
    # sencond the u-component
    ws.Copy(ws.casefull, ws.basename)
    ws.StringSet(ws.caseext, ".wind_u")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    # Print( casefull, 0 )
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
    # last the v-component
    ws.Copy(ws.casefull, ws.basename)
    ws.StringSet(ws.caseext, ".wind_v")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    # Print( casefull, 0 )
    # first: regrid to p_grid
    ws.GriddedFieldPRegrid(ws.gf3tmp, ws.p_grid, ws.gf3tmp, ws.interp_order, 1)
    # second: make 3D fields from 1D, then regrid to lat/lon_grid
    ws.GriddedFieldLatLonExpand(ws.gf3tmp, ws.gf3tmp)
    ws.GriddedFieldLatLonRegrid(
        ws.gf3tmp, ws.lat_grid, ws.lon_grid, ws.gf3tmp, ws.interp_order
    )
    # last, extract the data Tensor from the regridded GriddedField
    ws.FieldFromGriddedField(
        ws.wind_v_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
    )
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()


ws.forloop_agenda_solar = forloop_agenda_solar

ws.AgendaCreate("forloop_agenda_dust")


@arts_agenda
def forloop_agenda_dust(ws):
    # construct atmcase name III (Mars.LsXX.YY.dust-ZZ)
    ws.Extract(ws.casefull, ws.dustcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    # keep the casestring till dust and make upper-level folder name
    ws.Append(ws.basename, ws.atmcase)
    ws.StringSet(ws.caseext, "/")
    ws.Append(ws.basename, ws.caseext)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_solar)
    ws.nelemGet(ws.ncases, ws.solarcasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_dust = forloop_agenda_dust

ws.AgendaCreate("forloop_agenda_season")


@arts_agenda
def forloop_agenda_season(ws):
    # since we have different (number of) cases (Ne) depending on day/night, we
    # can't use a loop to process them, but need to set and call each separately.
    # hence we replace the timeloop (which would loops over 2 instances) by
    # explicit calls.
    # first: day
    # construct atmcase name I (Mars.LsXX)
    ws.Extract(ws.atmcase, ws.seasoncasearray, ws.forloop_index)
    # construct atmcase name II (Mars.LsXX.d/n)
    ws.Extract(ws.casefull, ws.timecasearray, 0)
    ws.Append(ws.atmcase, ws.casefull)
    # set to use the correct Ne case array
    ws.Copy(ws.necasearray, ws.daynearray)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_dust)
    ws.nelemGet(ws.ncases, ws.dustcasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
    # second: night
    # construct atmcase name I (Mars.LsXX)
    ws.Extract(ws.atmcase, ws.seasoncasearray, ws.forloop_index)
    # construct atmcase name II (Mars.LsXX.d/n)
    ws.Extract(ws.casefull, ws.timecasearray, 1)
    ws.Append(ws.atmcase, ws.casefull)
    # set to use the correct Ne case array
    ws.Copy(ws.necasearray, ws.nightnearray)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_dust)
    ws.nelemGet(ws.ncases, ws.dustcasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_season = forloop_agenda_season

ws.nelemGet(ws.ncases, ws.seasoncasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_season)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
