#
# Testing functionality (meeting format requirements, etc.) of supplemental
#  atmospheric scenario data.
#
# General test setup: reading in raw data (including a basic atmosphere),
#  expanding (where necessary), regridding (incl. extracting), executing
#  standard pre-RT calc internal test method atmfields_checkedCalc.
#
#
# This case is for Venus and specifically tests
#
#  - for the five Venus scenarions: Venus.spicav.night, Venus.spicav.night_cold,
#     Venus.vira.night, Venus.vira.day, Venus.vira.day_highlat
#  - electron densities (given: 1D cases, 2 per night scneario, 5 per day scenario)
#      - all cases in 1D with p_grid from altitude grid
#      - a single case expanded to 3D (assuming the other cases behave in the
#        same way) with p_grid from altitude grid
#  - wind (given: 1D cases with both horizontal components, 1 case for NS
#     component and 3 cases for EW component per scenario)
#      - all expanded to global 3D case with p_grid from altitude grid
#
# Jana Mendrok 2013-02-26

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_venus.arts")
# 1-dimensional atmosphere
ws.AtmosphereSet1D()
# We have to set the absorption lookup table interpolation order,
# since we are using wind fields. atmferlds_checkedCalc will otherwise throw an error.
ws.IndexSet(ws.abs_f_interp_order, 1)
ws.Tensor3Create("edensity_field")
ws.GriddedField3Create("gf3tmp")
ws.StringCreate("caseext")
ws.StringCreate("casefull")
ws.StringCreate("atmcase")
ws.IndexCreate("ncases")
ws.ArrayOfStringCreate("atmcasearray")
ws.ArrayOfStringCreate("necasearray")
ws.IndexCreate("interp_order")
ws.IndexSet(ws.interp_order, 1)
# set basic case folder
ws.StringCreate("basename")
ws.StringSet(ws.basename, "planets/Venus/MPS/")
# Array with case names
ws.ArrayOfStringCreate("nightcasearray")
ws.ArrayOfStringSet(
    ws.nightcasearray,
    ["Venus.spicav.night", "Venus.spicav.night_cold", "Venus.vira.night"],
)
ws.ArrayOfStringCreate("nightnearray")
ws.ArrayOfStringSet(ws.nightnearray, [".SZA.90-100.Ne", ".SZA.100-120.Ne"])
ws.ArrayOfStringCreate("daycasearray")
ws.ArrayOfStringSet(ws.daycasearray, ["Venus.vira.day", "Venus.vira.day_highlat"])
ws.ArrayOfStringCreate("daynearray")
ws.ArrayOfStringSet(
    ws.daynearray,
    [
        ".SZA.0-30.Ne",
        ".SZA.30-50.Ne",
        ".SZA.50-70.Ne",
        ".SZA.70-80.Ne",
        ".SZA.80-90.Ne",
    ],
)
ws.ArrayOfStringCreate("windcasearray")
ws.ArrayOfStringSet(ws.windcasearray, ["_max", "_mid", "_min"])
# the 3D geo grid to test
ws.VectorCreate("lat_grid3D")
ws.VectorCreate("lon_grid3D")
ws.VectorLinSpace(ws.lat_grid3D, -90.0, 90.0, 18.0)
ws.VectorLinSpace(ws.lon_grid3D, -20.0, 350.0, 18.0)
#####
# CASES A-C (night, 2 Ne profiles) and D-E (day, 5 Ne profiles)
#####
# we go with a foorloop through the different cases - night/day in an outer loop
# and the different Ne cases within in an inner loop.
ws.AgendaCreate("forloop_agenda_Ne")


@arts_agenda
def forloop_agenda_Ne(ws):
    # construct the full case name
    ws.Copy(ws.casefull, ws.atmcase)
    ws.Extract(ws.caseext, ws.necasearray, ws.forloop_index)
    ws.Append(ws.casefull, ws.caseext)
    ws.Print(ws.casefull, 0)
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

ws.AgendaCreate("forloop_agenda_wind")


@arts_agenda
def forloop_agenda_wind(ws):
    # construct the full case name
    ws.Copy(ws.atmcase, ws.casefull)
    ws.Extract(ws.caseext, ws.windcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.caseext)
    ws.Print(ws.atmcase, 0)
    # readin in raw data
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
        ws.wind_v_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf3tmp
    )
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()


ws.forloop_agenda_wind = forloop_agenda_wind

ws.AgendaCreate("forloop_agenda_scen")


@arts_agenda
def forloop_agenda_scen(ws):
    # construct atmscen name
    ws.Copy(ws.atmcase, ws.basename)
    ws.Extract(ws.casefull, ws.atmcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    ws.StringSet(ws.caseext, "/")
    ws.Append(ws.atmcase, ws.caseext)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Print(ws.atmcase, 0)
    # we manually select a minumim set of basic atm data (main atm constituents)
    ws.abs_speciesSet(species=["CO2"])
    ws.AtmRawRead(basename=ws.atmcase)
    ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
    ws.AtmFieldsCalc(vmr_zeropadding=1)
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
    ws.Copy(ws.casefull, ws.atmcase)
    ws.Extract(ws.caseext, ws.necasearray, 0)
    ws.Append(ws.casefull, ws.caseext)
    ws.Print(ws.casefull, 0)
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
    # first the v-component (as this is not changing within one scenario)
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".wind_v")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.Print(ws.casefull, 0)
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
    # now the u-component. for the 3 cases per scenario, we again use a forloop.
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".wind_u")
    ws.Append(ws.casefull, ws.caseext)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_wind)
    ws.nelemGet(ws.ncases, ws.windcasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_scen = forloop_agenda_scen

ws.Copy(ws.forloop_agenda, ws.forloop_agenda_scen)
ws.Copy(ws.atmcasearray, ws.nightcasearray)
ws.Copy(ws.necasearray, ws.nightnearray)
ws.nelemGet(ws.ncases, ws.atmcasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
ws.Copy(ws.atmcasearray, ws.daycasearray)
ws.Copy(ws.necasearray, ws.daynearray)
ws.nelemGet(ws.ncases, ws.atmcasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
