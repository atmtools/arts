#
# Testing functionality (meeting format requirements, etc.) of basic atmospheric
#  scenario data.
#
# General test setup: reading in raw data, regridding to common p-grid (1D),
#  executing standard pre-RT calc internal test method atmfields_checkedCalc.
#
#
# This case is for Mars and specifically tests
#
# (CASES 1-72)
#  - 72 Mars scenarions: 4seasons x 2daytimes x 3dustloads x 3solaractivities
#  - t, z, and all abs species vmr (the ones that follow the basename convention
#     are caught by abs_speciesDefineAllInScenario and AtmRawRead; others are
#     derived manually) in the case folder (no abs species in "unused/"
#     subfolder).
#  - regridding to a pressure grid taken from the read-in altitude grid
#
# Jana Mendrok 2013-02-26

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_mars.arts")
# 1-dimensional atmosphere
ws.AtmosphereSet1D()
ws.GriddedField3Create("gf3tmp")
ws.StringCreate("caseext")
ws.StringCreate("casefull")
ws.StringCreate("atmcase")
ws.IndexCreate("ncases")
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
# we go with several nested foorloop through the different cases.
#  All those cases have identical abs species to process.
#  Order of agenda definitions has to be inverse from their execution (as we can
#  only copy an agenda AFTER we have defined it).
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
    ws.Print(ws.basename, 0)
    # derive absspecies with standard name from scenario
    ws.abs_speciesDefineAllInScenario(basename=ws.basename)
    # WriteXMLIndexed( "ascii", forloop_index,
    #                 abs_species, "TestAtmScen_Mars_allInScen.abs_species" )
    ws.AtmRawRead(basename=ws.basename)
    # adding species or variants that do not follow the general naming convention
    ws.abs_speciesAdd(species=["CH4"])
    ws.Copy(ws.casefull, ws.basename)
    ws.StringSet(ws.caseext, ".CH4_high")
    ws.Append(ws.casefull, ws.caseext)
    # Print( casefull, 0 )
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.Append(ws.vmr_field_raw, ws.gf3tmp)
    ws.abs_speciesAdd(species=["H2O-162"])
    ws.Copy(ws.casefull, ws.basename)
    ws.StringSet(ws.caseext, ".H2O-162")
    ws.Append(ws.casefull, ws.caseext)
    # Print( casefull, 0 )
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.Append(ws.vmr_field_raw, ws.gf3tmp)
    # now derive common p_grid and regrid atm fields to this
    ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
    ws.AtmFieldsCalc()
    ws.Extract(ws.z_surface, ws.z_field, 0)
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()
    # WriteXML( "ascii", p_grid )
    # WriteXML( "ascii", z_field )
    # WriteXML( "ascii", t_field )
    # WriteXML( "ascii", vmr_field_raw )
    ws.WriteXMLIndexed("ascii", ws.forloop_index, ws.vmr_field)


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

ws.AgendaCreate("forloop_agenda_time")


@arts_agenda
def forloop_agenda_time(ws):
    # construct atmcase name II (Mars.LsXX.d/n)
    ws.Extract(ws.casefull, ws.timecasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_dust)
    ws.nelemGet(ws.ncases, ws.dustcasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_time = forloop_agenda_time

ws.AgendaCreate("forloop_agenda_season")


@arts_agenda
def forloop_agenda_season(ws):
    # construct atmcase name I (Mars.LsXX)
    ws.Extract(ws.casefull, ws.seasoncasearray, ws.forloop_index)
    ws.Copy(ws.atmcase, ws.casefull)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_time)
    ws.nelemGet(ws.ncases, ws.timecasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_season = forloop_agenda_season

ws.nelemGet(ws.ncases, ws.seasoncasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_season)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
